import re

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.shortcuts import render, redirect, get_object_or_404
from django.db.models import Q
from django.http import HttpResponse, Http404, HttpResponseRedirect
from django.conf import settings

import os

from django.views.generic import TemplateView, ListView, DetailView

from main.models import Compound, Plant
from core.utils.QueryHandler import query_to_df, df_to_sdf, update_sdf


class AboutView(TemplateView):
    extra_context = {
        'n_compounds': Compound.objects.all().count(),
        'n_plants': Plant.objects.all().count()
    }
    template_name = 'main/about.html'


class QueryResultListView(ListView):
    paginate_by = 10
    template_name = 'main/results.html'
    context_object_name = 'compounds'

    def get_queryset(self):
        query = self.request.GET.get('q')
        if query:
            return Compound.objects.filter(Q(id=query)
                                           | Q(PID=query)
                                           | Q(Smiles=query)
                                           | Q(Molecular_Formula=query)
                                           | Q(plants__name__iexact=query)).distinct()
        else:
            return Compound.objects.none()


class PlantCompoundsListView(ListView):
    template_name = 'main/plant.html'
    context_object_name = 'compounds'

    def get_queryset(self):
        self.plant = get_object_or_404(Plant, id=self.kwargs['pk'])
        return Compound.objects.filter(plants=self.plant).order_by('id')

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['plant'] = self.plant
        return context


class CompoundDetailView(DetailView):
    model = Compound
    template_name = 'main/compound.html'
    context_object_name = 'compound'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        lipinski = (self.get_object().H_Bond_Donors > 5) + (self.get_object().H_Bond_Acceptors > 10) + \
                   (self.get_object().Molecular_Weight >= 500) + (self.get_object().logP > 5)
        violation_color = ['#2ECC71', '#ABEBC6', '#FADBD8', '#E74C3C', '#B03A2E']
        context['lipinski'] = lipinski
        context['violation_color'] = violation_color
        return context


# view for download a query results
def download_file(request, search):
    # make search term as is in the PID column of database
    if search.isnumeric():
        search = 'Phytochem_' + search.zfill(6)
    # temporary path for storing download files
    temppath = os.path.join(settings.MEDIA_ROOT, 'tempDownloadFiles/')

    if os.path.isdir(temppath):
        oldfiles = os.listdir(temppath)
        for file in oldfiles:
            try:
                os.remove(temppath + file)
            except FileNotFoundError:
                pass
    else:
        os.makedirs(temppath, exist_ok=True)

    # search in the database columns PID, Smiles and
    # Molecular_Formula
    compounds = Compound.objects.filter(
        Q(PID=search) | Q(Smiles=search) | Q(Molecular_Formula=search) | Q(plants__name__iexact=search))
    compounds_df = query_to_df(compounds)
    dw_file = os.path.join(temppath, search + '.sdf')
    df_to_sdf(compounds_df, dw_file)
    return prepare_download(dw_file)


@login_required
# View for downloading all file
def download_all_file(request):
    user_email = request.user.email
    tld = re.search("@\w+\.([\w]+)([\.\w]*)", user_email).group(1)
    name = 'all_data.sdf'
    file_path = os.path.join(settings.MEDIA_ROOT, name)
    if tld == 'ac' or tld == 'edu':
        return prepare_download(file_path)
    else:
        messages.warning(request, 'You must login with your institutional mail to download the dataset')
        return HttpResponseRedirect(request.META.get('HTTP_REFERER'))


# not view
def prepare_download(path):
    if not os.path.exists(path):
        update_sdf()
    if os.path.exists(path):
        with open(path, 'rb') as fh:
            response = HttpResponse(fh.read(), content_type='application/octet-stream')
            response['Content-Disposition'] = 'inline; filename=' + os.path.basename(path)
            return response
    raise Http404

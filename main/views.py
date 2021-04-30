from django.contrib import messages
from django.contrib.auth.mixins import LoginRequiredMixin
from django.shortcuts import get_object_or_404
from django.db.models import Q
from django.http import HttpResponse, HttpResponseForbidden, FileResponse
from django.conf import settings
from django.views import View
from django.views.generic import TemplateView, ListView, DetailView

import re
from io import StringIO
import pandas as pd
from rdkit.Chem import PandasTools

from main.models import Compound, Plant


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


class SDFResponse(HttpResponse):
    def __init__(self, data, output_name):
        output = StringIO()
        compounds_df = pd.DataFrame(list(data.values())).drop('id', axis=1)
        PandasTools.AddMoleculeColumnToFrame(compounds_df, 'Smiles', 'ROMol', includeFingerprints=True)
        PandasTools.WriteSDF(compounds_df, output, molColName='ROMol', idName='PID',
                             properties=list(compounds_df.columns))

        mimetype = 'text/plain'
        file_ext = 'sdf'
        output.seek(0)
        super(SDFResponse, self).__init__(content=output.getvalue(),
                                          content_type=mimetype)
        self['Content-Disposition'] = 'attachment;filename="%s.%s"' % \
                                      (output_name.replace('"', '\"'), file_ext)


class FileDownloadView(View):
    def get(self, request, search):
        compounds = Compound.objects.filter(Q(id=search)
                                            | Q(PID=search)
                                            | Q(Smiles=search)
                                            | Q(Molecular_Formula=search)
                                            | Q(plants__name__iexact=search))
        return SDFResponse(compounds, search)


class FullDownloadView(LoginRequiredMixin, View):
    def get(self, request):
        academic_mail = re.search(r"@\w+\.([\w]+)([\.\w]*)", request.user.email).group(1) in ['ac', 'edu']
        if not request.user.is_superuser and not academic_mail:
            messages.warning(request, 'You must login with your institutional mail to download the dataset')
            return HttpResponseForbidden()
        absolute_path = '{}/{}'.format(settings.MEDIA_ROOT, 'all_data.sdf')
        response = FileResponse(open(absolute_path, 'rb'), as_attachment=True, content_type='text/plain')
        return response

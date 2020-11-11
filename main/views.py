from django.shortcuts import render, redirect
from django.db.models import Q
from django.http import FileResponse, HttpResponse, Http404
from django.conf import settings

import os
from main.models import Compound, Plant
from utils.QueryHandler import query_to_df, df_to_sdf, df_to_pdb, df_to_mol


def index(request):
    context = {
        'n_elements': Compound.objects.all().count()
    }
    return render(request, 'main/index.html', context=context)


def results(request):
    search = request.GET['search']
    if len(search) != 0:
        if search.isnumeric():
            search = 'Phytochem_' + search.zfill(6)
        compounds = Compound.objects.filter(Q(PID=search) | Q(Smiles=search) | Q(Molecular_Formula=search))
        # print(search, compounds)
        context = {
            'title': 'Search results',
            'compounds': compounds
        }
        return render(request, 'main/results.html', context=context)
    return redirect('/')


def plant(request, id):
    plant = Plant.objects.get(id=id)
    compounds = plant.compound_set.all()
    context = {
        'title': 'Plant | Details',
        'plant': plant,
        'compounds': compounds
    }
    return render(request, 'main/plant.html', context=context)


def compound(request, id):
    compound = Compound.objects.get(id=id)
    print(Compound)
    plants = compound.plants.all()
    context = {
        'compound': compound,
        'plants': plants
    }
    return render(request, 'main/compound.html', context=context)


# view for download a query results
def download_file(request):
    search = request.GET['search']
    filetype = request.GET['filetype']
    # make search term as is in the PID column of database
    if search.isnumeric():
        search = 'Phytochem_' + search.zfill(6)

    # temporary path for storing download files
    temppath = os.path.join(settings.MEDIA_ROOT, 'tempDownloadFiles/')

    try:
        os.removedirs(temppath)
        os.makedirs(temppath, exist_ok=True)
    except FileNotFoundError:
        os.makedirs(temppath, exist_ok=True)

    # search in the database columns PID, Smiles and
    # Molecular_Formula
    compounds = Compound.objects.filter(Q(PID=search) |
                                        Q(Smiles=search) |
                                        Q(Molecular_Formula=search))
    compounds_df = query_to_df(compounds)

    if filetype == 'sdf':
        dw_file = os.path.join(temppath, search + '.sdf')
        df_to_sdf(compounds_df, dw_file)
    elif filetype == 'pdb':
        dw_file = os.path.join(temppath, search + '.pdb')
        df_to_pdb(compounds_df, dw_file)
    elif filetype == 'mol':
        dw_file = os.path.join(temppath, search + '.mol')
        df_to_mol(compounds_df, dw_file)
    return prepare_download(dw_file)


# View for downloading all file
def download_all_file(request):
    name = 'all_data.sdf'
    file_path = os.path.join(settings.MEDIA_ROOT, name)
    return prepare_download(file_path)


# not view
def prepare_download(path):
    if os.path.exists(path):
        with open(path, 'rb') as fh:
            response = HttpResponse(fh.read(), content_type='application/octet-stream')
            response['Content-Disposition'] = 'inline; filename=' + os.path.basename(path)
            return response
    raise Http404

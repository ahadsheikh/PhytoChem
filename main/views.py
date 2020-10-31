from django.shortcuts import render
from django.db.models import Q
from django.http import FileResponse

import os
import tempfile
from main.models import Compound
from utils.QueryHandler import query_to_df, df_to_sdf, df_to_pdb, df_to_mol


def index(request):
    context = {
        'title': "Home"
    }
    return render(request, 'main/index.html', context=context)


def results(request):
    search = request.GET['search']
    if search.isnumeric():
        search = 'Phytochem_' + search.zfill(6)
    compounds = Compound.objects.filter(Q(PID=search) | Q(Smiles=search) | Q(Molecular_Formula=search))
    print(search, compounds)
    context = {
        'title': 'Search results',
        'compounds': compounds
    }
    return render(request, 'main/results.html', context=context)


def download_file(request):
    try:
        search = request.GET['search']
        filetype = request.GET['filetype']
        # make search term as is in the PID column of database
        if search.isnumeric():
            search = 'Phytochem_' + search.zfill(6)
        # search in the database columns PID, Smiles and
        # Molecular_Formula
        compounds = Compound.objects.filter(Q(PID=search) |
                                            Q(Smiles=search) |
                                            Q(Molecular_Formula=search))
        compounds_df = query_to_df(compounds)

        # create a file based on file extension
        tmp = tempfile.NamedTemporaryFile(suffix=".{}".format(filetype),
                                          prefix="pc_{}_".format(search),
                                          delete=False)
        if filetype == 'sdf':
            df_to_sdf(compounds_df, tmp.name)
        elif filetype == 'pdb':
            df_to_pdb(compounds_df, tmp.name)
        elif filetype == 'mol':
            df_to_mol(compounds_df, tmp.name)
        response = FileResponse(open(tmp.name, 'rb'))
        return response
    finally:
        os.remove(tmp.name)

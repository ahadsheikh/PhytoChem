from django.shortcuts import render
from django.db.models import Q
from django.urls import reverse

from main.models import Compound
import os
import tempfile

# necessary imports for creating files
from django.http import FileResponse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools

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
        if search.isnumeric():
            search = 'Phytochem_' + search.zfill(6)
        compounds = Compound.objects.filter(Q(PID=search) | Q(Smiles=search) | Q(Molecular_Formula=search))
        compounds_df = query_to_df(compounds)

        if filetype == 'sdf':
            tmp = tempfile.NamedTemporaryFile(suffix=".sdf", prefix="pc_{}_".format(search), delete=False)
            df_to_sdf(compounds_df, tmp)
        elif filetype == 'pdb':
            tmp = tempfile.NamedTemporaryFile(suffix=".pdb", prefix="pc_{}_".format(search), delete=False)
            df_to_pdb(compounds_df, tmp)
        elif filetype == 'mol':
            tmp = tempfile.NamedTemporaryFile(suffix=".mol", prefix="pc_{}_".format(search), delete=False)
            df_to_mol(compounds_df, tmp)

        response = FileResponse(open(tmp.name, 'rb'))
        return response
    finally:
        os.remove(tmp.name)

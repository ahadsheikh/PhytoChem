from django.shortcuts import render
from django.db.models import Q
from django.urls import reverse

from main.models import Compound
import os
import tempfile
from django.http import FileResponse

# necessary imports for creating files
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import Crippen
import glob
import pandas as pd


def index(request):
    context = {
        'title': "Home"
    }
    return render(request, 'main/index.html', context=context)


def results(request):
    search = request.GET['search']
    if search.isnumeric():
        search = 'Phytochem_' + search.zfill(5)
    compounds = Compound.objects.filter(Q(PID=search) | Q(Smiles=search) | Q(Molecular_Formula=search))
    context = {
        'title': 'Search results',
        'compounds': compounds
    }
    return render(request, 'main/results.html', context=context)


def download_file(request):
    try:
        pname = "pcs"
        search = request.GET['search']
        filetype = request.GET['filetype']
        if search.isnumeric():
            search = 'Phytochem_' + search.zfill(5)
        compounds = Compound.objects.filter(Q(PID=search) | Q(Smiles=search) | Q(Molecular_Formula=search))
        # Dataframe to write calculations of each compounds
        compounds_df = pd.DataFrame(
            columns=['ID', 'Smiles', 'Molecular_Formula', 'Molecular_Weight', 'H_Bond_Acceptors',
                     'H_Bond_Donors', 'Molar_Refractivity', 'TPSA'])
        for compound in compounds:
            compounds_df = compounds_df.append({
                'ID': compound.PID,
                'Smiles': compound.Smiles,
                'Molecular_Formula': Compound.Molecular_Formula,
                'Molecular_Weight': Compound.Molecular_Weight,
                'H_Bond_Acceptors': Compound.H_Bond_Acceptors,
                'H_Bond_Donors': Compound.H_Bond_Donors,
                'Molar_Refractivity': Compound.Molar_Refractivity,
                'TPSA': Compound.TPSA,
            }, ignore_index=True)
        pname = compounds[0].Molecular_Formula
        PandasTools.AddMoleculeColumnToFrame(compounds_df, 'Smiles', 'ROMol', includeFingerprints=True)
        if filetype == 'sdf':
            tmp = tempfile.NamedTemporaryFile(suffix=".sdf", prefix="pc_{}_".format(pname), delete=False)
            with open(tmp.name, 'w') as fi:
                PandasTools.WriteSDF(compounds_df, fi, molColName='ROMol', idName='ID',
                                     properties=list(compounds_df.columns))
        elif filetype == 'pdb':
            tmp = tempfile.NamedTemporaryFile(suffix=".pdb", prefix="pc_{}_".format(pname), delete=False)
            with open(tmp.name, 'w') as fi:
                pdbwriter = Chem.PDBWriter(fi)
                for compound in compounds:
                    mol = Chem.MolFromSmiles(compound.Smiles)
                    pdbwriter.write(mol)
                pdbwriter.close()
        elif filetype == 'mol':
            tmp = tempfile.NamedTemporaryFile(suffix=".mol", prefix="pc_{}_".format(pname), delete=False)
            with open(tmp.name, 'w+') as fi:
                for compound in compounds:
                    mol = Chem.MolFromSmiles(compound.Smiles)
                    molblock = Chem.MolToMolBlock(mol)
                    fi.write(molblock)
        response = FileResponse(open(tmp.name, 'rb'))
        return response
    finally:
        os.remove(tmp.name)

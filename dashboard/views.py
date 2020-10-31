from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
from django.shortcuts import render, redirect
from django.contrib import messages

import os
import pandas as pd

from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import Crippen

from dashboard.forms import UploadFileForm
from main.models import Compound
from utils.QueryHandler import update_sdf


@login_required(redirect_field_name='next')
def dash_index(request):
    context = {
        'title': 'Dashboard'
    }
    return render(request, 'dashboard/dash.html', context=context)


@login_required(redirect_field_name='next')
def upload(request):
    if request.method == 'POST':
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            handle_file(request.FILES['file'])
            messages.success(request, "Upload File Successfully")
            return redirect('dashboard')

    messages.success(request, "Upload File Not Successfully")
    return HttpResponse(
        "<h2>You are not permitted.</h2>"
    )


# Not path view
def handle_file(f):
    if not os.path.isdir('media/upload'):
        os.mkdir('media')
        os.mkdir('media/upload')
    with open('media/upload/' + f.__str__(), 'wb') as des:
        for chunk in f.chunks():
            des.write(chunk)

    sdf_file = 'media/upload/' + f.__str__()
    sdf = Chem.SDMolSupplier(sdf_file)  # read sdf
    compounds_df = pd.DataFrame(
        columns=['Smiles', 'Molecular_Formula', 'Molecular_Weight', 'H_Bond_Acceptors',
                 'H_Bond_Donors', 'Molar_Refractivity', 'TPSA'])
    for mol in sdf:
        smiles = Chem.MolToSmiles(mol)  # get smiles
        if not Compound.objects.filter(Smiles=smiles):
            molecular_formula = Chem.rdMolDescriptors.CalcMolFormula(mol)  # formula
            molecular_weight = Chem.rdMolDescriptors.CalcExactMolWt(mol)  # weight
            hba = Chem.rdMolDescriptors.CalcNumHBA(mol)  # h bond acceptor
            hbd = Chem.rdMolDescriptors.CalcNumHBD(mol)  # h bond donor
            molar_refractivity = Chem.Crippen.MolMR(mol)  # molar refractivity
            tpsa = Chem.rdMolDescriptors.CalcTPSA(mol)  # tpsa
            compounds_df = compounds_df.append({  # write this row to dataframe
                'Smiles': smiles,
                'Molecular_Formula': molecular_formula,
                'Molecular_Weight': molecular_weight,
                'H_Bond_Acceptors': hba,
                'H_Bond_Donors': hbd,
                'Molar_Refractivity': molar_refractivity,
                'TPSA': tpsa
            }, ignore_index=True)
    PandasTools.AddMoleculeColumnToFrame(compounds_df, 'Smiles', 'ROMol', includeFingerprints=True)
    compounds_df.drop_duplicates(subset="Smiles", keep=False, inplace=True)  # drop duplicate by smiles
    compound_len = Compound.objects.all().count()  # get current compounds in database
    counter = 0
    for i, compound in compounds_df.iterrows():
        counter += 1
        Compound.objects.create(
            PID='Phytochem_' + str(compound_len + counter).zfill(6),
            Smiles=compound.Smiles,
            Molecular_Formula=compound.Molecular_Formula,
            Molecular_Weight=compound.Molecular_Weight,
            H_Bond_Acceptors=compound.H_Bond_Acceptors,
            H_Bond_Donors=compound.H_Bond_Donors,
            Molar_Refractivity=compound.Molar_Refractivity,
            TPSA=compound.TPSA,
            ROMol=compound.ROMol
        )
    update_sdf()

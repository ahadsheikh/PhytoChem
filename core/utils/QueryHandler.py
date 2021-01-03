import os

import pandas as pd
from django.core.exceptions import ValidationError
from rdkit import Chem
from rdkit.Chem import PandasTools, Crippen
from main.models import Compound, Plant
from bs4 import BeautifulSoup


def query_to_df(queryset):
    # Dataframe to write calculations of each compounds
    compounds_df = pd.DataFrame(list(queryset.values())).drop('id', axis=1)
    PandasTools.AddMoleculeColumnToFrame(compounds_df, 'Smiles', 'ROMol', includeFingerprints=True)
    return compounds_df


def df_to_sdf(compounds_df, file):
    with open(file, 'w') as fi:
        PandasTools.WriteSDF(compounds_df, fi, molColName='ROMol', idName='PID',
                             properties=list(compounds_df.columns))


def get_src_from_image_tag(html):
    soup = BeautifulSoup(str(html), "html.parser")
    return soup.img['src']


def update_sdf():
    compounds_df = pd.DataFrame(list(Compound.objects.all().values())).drop(['id', 'created_at', 'updated_at'], axis=1)
    PandasTools.AddMoleculeColumnToFrame(compounds_df, 'Smiles', 'ROMol', includeFingerprints=True)
    if not os.path.exists('media'):
        os.makedirs('media')
    df_to_sdf(compounds_df, 'media/all_data.sdf')


def update_db_from_df(compounds_df, plant=None):
    compound_len = Compound.objects.all().count()  # get current compounds in database
    counter = 0

    for i, compound in compounds_df.iterrows():
        if not Compound.objects.filter(Smiles=compound.Smiles):
            counter += 1
            cur_compound = Compound.objects.create(
                PID='Phytochem_' + str(compound_len + counter).zfill(6),
                Smiles=compound.Smiles,
                Molecular_Formula=compound.Molecular_Formula,
                Molecular_Weight=compound.Molecular_Weight,
                H_Bond_Acceptors=compound.H_Bond_Acceptors,
                H_Bond_Donors=compound.H_Bond_Donors,
                Molar_Refractivity=compound.Molar_Refractivity,
                TPSA=compound.TPSA,
                logP=compound.logP,
                ROMol=compound.ROMol
            )
        else:
            cur_compound = Compound.objects.get(Smiles=compound.Smiles)
        if plant:
            plant, _ = Plant.objects.get_or_create(name=plant)
            cur_compound.plants.add(plant)
    update_sdf()


def handle_new_sdf(path, plant=None, change_db=True):
    sdf = Chem.SDMolSupplier(path)  # read sdf
    compounds_df = pd.DataFrame(
        columns=['Smiles', 'Molecular_Formula', 'Molecular_Weight', 'H_Bond_Acceptors',
                 'H_Bond_Donors', 'Molar_Refractivity', 'TPSA', 'logP'])
    for mol in sdf:
        smiles = Chem.MolToSmiles(mol)  # get smiles
        molecular_formula = Chem.rdMolDescriptors.CalcMolFormula(mol)  # formula
        molecular_weight = Chem.rdMolDescriptors.CalcExactMolWt(mol)  # weight
        hba = Chem.rdMolDescriptors.CalcNumHBA(mol)  # h bond acceptor
        hbd = Chem.rdMolDescriptors.CalcNumHBD(mol)  # h bond donor
        molar_refractivity = Chem.Crippen.MolMR(mol)  # molar refractivity
        tpsa = Chem.rdMolDescriptors.CalcTPSA(mol)  # tpsa
        logp = Chem.Crippen.MolLogP(mol)
        compounds_df = compounds_df.append({  # write this row to dataframe
            'Smiles': smiles,
            'Molecular_Formula': molecular_formula,
            'Molecular_Weight': molecular_weight,
            'H_Bond_Acceptors': hba,
            'H_Bond_Donors': hbd,
            'Molar_Refractivity': molar_refractivity,
            'TPSA': tpsa,
            'logP': logp
        }, ignore_index=True)
    PandasTools.AddMoleculeColumnToFrame(compounds_df, 'Smiles', 'ROMol', includeFingerprints=True)
    compounds_df.ROMol = compounds_df.ROMol.apply(get_src_from_image_tag)
    if change_db:
        update_db_from_df(compounds_df, plant)
        update_sdf()
    else:
        return compounds_df


def validate_sdf(path):
    df = PandasTools.LoadSDF(path)
    if df.empty:
        raise ValidationError(
            'File Validation Failed',
            params={'file': path},
        )

import csv
from io import StringIO
import pandas as pd
import requests
from rdkit import Chem
import glob
import os
import re
from bs4 import BeautifulSoup
from rdkit.Chem import PandasTools
from main.models import Compound, Plant


def update_sdf():
    compounds_df = pd.DataFrame(list(Compound.objects.all().values()))
    if not compounds_df.isnull:
        compounds_df = compounds_df.drop(['id', 'created_at', 'updated_at'], axis=1)
    PandasTools.AddMoleculeColumnToFrame(compounds_df, 'Smiles', 'ROMol', includeFingerprints=True)
    if not os.path.exists('media'):
        os.makedirs('media')
    with open('media/all_data.sdf', 'w') as fi:
        PandasTools.WriteSDF(compounds_df, fi, molColName='ROMol', idName='PID',
                             properties=list(compounds_df.columns))


def single_smiles_scrape(smiles: str) -> str:
    # Give SwissADME data from  a smiles
    url = "http://www.swissadme.ch/index.php"
    s = requests.post(url, data={'smiles': smiles})

    soup = BeautifulSoup(s.content, 'html.parser')
    csv_div = soup.find('div', attrs={
        "style": "float: right; width: 250px; height: 30px; text-align: left; font: 13pt  Helvetica,Verdana, "
                 "sans-serif;"})
    csv_link = "http://www.swissadme.ch/" + csv_div.a['href']
    csv_data = requests.get(csv_link)
    return csv_data.text


def swisscsv_to_db(cur_compound_no, smiles, mol):
    csvdata = StringIO(single_smiles_scrape(smiles))
    csv_reader = csv.reader(csvdata, delimiter=',')

    lc = 0
    for row in csv_reader:
        if lc == 1:
            compound = Compound()
            compound.PID = 'Phytochem_{}'.format(str(cur_compound_no).zfill(6))
            print("adding", compound.PID)
            compound.Smiles = smiles
            compound.ROMol = get_src_from_image_tag(str(mol))  # image of the molecule

            if row[0].startswith('Molecule'):
                row[0] = None
            compound.Molecule_Name = row[0]
            # Physicochemical Properties
            compound.Molecular_Formula = row[2]
            compound.Molecular_Weight = row[3]
            compound.Heavy_Atoms = row[4]
            compound.Arom_Heavy_Atoms = row[5]
            compound.Fraction_Csp3 = row[6]
            compound.Rotatable_bonds = row[7]
            compound.H_Bond_Acceptors = row[8]
            compound.H_Bond_Donors = row[9]
            compound.Molar_Refractivity = row[10]
            compound.TPSA = row[11]

            # Lipophilicity
            compound.iLOGP = row[12]

            # Pharmacokinetics
            compound.GI_absorption = row[30]
            compound.BBB_permeant = row[31]
            compound.P_gp_substrate = row[32]
            compound.CYP1A2_inhibitor = row[33]
            compound.CYP2C19_inhibitor = row[34]
            compound.CYP2C9_inhibitor = row[35]
            compound.CYP2D6_inhibitor = row[36]
            compound.CYP3A4_inhibitor = row[37]
            compound.LogKp = row[38]

            # Druglikeness
            compound.Lipinski = row[39]
            compound.Ghose = row[40]
            compound.Veber = row[41]
            compound.Egan = row[42]
            compound.Bioavailability_Score = row[44]

            # Medicinal Chemistry
            compound.PAINS = row[45]
            compound.Brenk = row[46]
            compound.Leadlikeness = row[47]
            compound.Synthetic_accessibility = row[48]

            compound.save()

            return compound

        lc += 1


def get_src_from_image_tag(html):
    soup = BeautifulSoup(html, "html.parser")
    return soup.img['src']


def save_sdf(plant, path):
    sdf = Chem.SDMolSupplier(path)  # read sdf
    cur_compound_no = Compound.objects.order_by('PID').last()
    for mol in sdf:
        smiles = Chem.MolToSmiles(mol)  # get smiles
        # check if this compound already exists
        cur_compound = Compound.objects.filter(Smiles=smiles)
        if not cur_compound:
            # save the values in database
            cur_compound_no += 1
            cur_compound = [swisscsv_to_db(cur_compound_no, smiles, mol)]
        cur_compound = cur_compound[0]
        try:
            # if this plant already exists, just append to the plants_object list
            plant = Plant.objects.get(name=plant)
        except Plant.DoesNotExist:
            plant = Plant.objects.create(name=plant)
        cur_compound.plants.add(plant)


# option 0: single sdf
# option 1: plant name / sdf structured data
# option 2: plant name / *** / sdf structured data
def save_to_db(plant, path, filetype):
    if filetype == 0:
        save_sdf(plant, path)
    else:
        for file in glob.iglob('{}/**//*.sdf'.format(path), recursive=True):
            try:
                parent = os.path.dirname(file).split(os.sep)[-filetype]
                plant = re.search('([A-Za-z]+( [A-Za-z]+)?)', parent).group(0)
                save_sdf(plant, file)
            except Exception as e:
                print(e)
    update_sdf()

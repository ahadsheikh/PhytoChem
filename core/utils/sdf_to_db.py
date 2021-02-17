from rdkit import Chem
from rdkit.Chem import Crippen
import glob
import os
import re
from bs4 import BeautifulSoup

from core.utils.QueryHandler import update_sdf
from main.models import Compound, Plant


def dict_to_db(comp_dict):
    data = comp_dict
    data_len = Compound.objects.all().count()
    i = 0
    for k, v in data.items():
        # get data from the dictionary
        smiles = k
        molecular_formula = v['molecular_formula']
        molecular_weight = v['molecular_weight']
        h_bond_acceptors = v['hba']
        h_bond_donors = v['hbd']
        molar_refractivity = v['molar_refractivity']
        tpsa = v['tpsa']
        logp = v['logp']
        romol = v['romol']
        plants = v.get('plant', {})

        # check if this compound already exists
        cur_compound = Compound.objects.filter(Smiles=smiles)
        if not cur_compound:
            # save the values in database
            compound = Compound()
            i += 1
            compound.PID = 'Phytochem_{}'.format(str(data_len + i).zfill(6))
            compound.Smiles = smiles
            compound.Molecular_Formula = molecular_formula
            compound.Molecular_Weight = molecular_weight
            compound.H_Bond_Acceptors = h_bond_acceptors
            compound.H_Bond_Donors = h_bond_donors
            compound.Molar_Refractivity = molar_refractivity
            compound.TPSA = tpsa
            compound.logP = logp
            compound.ROMol = romol
            compound.save()
            cur_compound = [compound]
        cur_compound = cur_compound[0]
        for plant in plants:
            try:
                # if this plant already exists, just append to the plants_object list
                plant = Plant.objects.get(name=plant)
            except Plant.DoesNotExist:
                plant = Plant.objects.create(name=plant)
            cur_compound.plants.add(plant)

        if i % 100 == 0 and i > 0:
            print('{} data processed'.format(i))
    update_sdf()


def get_src_from_image_tag(html):
    soup = BeautifulSoup(html, "html.parser")
    return soup.img['src']


# option 0: single sdf
# option 1: plant name / sdf structured data
# option 2: plant name / *** / sdf structured data
def save_sdf(path, option):
    smiles_dict = dict()
    i = 0
    for file in glob.iglob('{}/**//*.sdf'.format(path), recursive=True):
        sdf = Chem.SDMolSupplier(file)  # read sdf
        for mol in sdf:
            smiles = Chem.MolToSmiles(mol)  # get smiles
            if smiles_dict.get(smiles, -1) == -1:
                i += 1
                smiles_dict[smiles] = {
                    'molecular_formula': Chem.rdMolDescriptors.CalcMolFormula(mol),  # formula
                    'molecular_weight': Chem.rdMolDescriptors.CalcExactMolWt(mol),  # weight
                    'hba': Chem.rdMolDescriptors.CalcNumHBA(mol),  # h bond acceptor
                    'hbd': Chem.rdMolDescriptors.CalcNumHBD(mol),  # h bond donor
                    'molar_refractivity': Chem.Crippen.MolMR(mol),  # molar refractivity
                    'tpsa': Chem.rdMolDescriptors.CalcTPSA(mol),  # tpsa
                    'logp': Chem.Crippen.MolLogP(mol),  # logP
                    'romol': get_src_from_image_tag(str(mol))  # image of the molecule
                }
                if option:
                    parent = os.path.dirname(file).split(os.sep)[-option]
                    plant = re.search('([A-Za-z]+( [A-Za-z]+)?)', parent).group(0)
                    smiles_dict[smiles]['plant'] = [plant]
            else:
                if option:
                    parent = os.path.dirname(file).split(os.sep)[-option]
                    plant = re.search('([A-Za-z]+( [A-Za-z]+)?)', parent).group(0)
                    if plant not in smiles_dict[smiles]['plant']:
                        smiles_dict[smiles]['plant'].append(plant)
            if i % 100 == 0 and i > 0:
                print('{} compounds processed'.format(i))
    dict_to_db(smiles_dict)

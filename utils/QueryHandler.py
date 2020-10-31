import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
from main.models import Compound


def query_to_df(queryset):
    # Dataframe to write calculations of each compounds
    compounds_df = pd.DataFrame(list(queryset.values())).drop('id', axis=1)
    PandasTools.AddMoleculeColumnToFrame(compounds_df, 'Smiles', 'ROMol', includeFingerprints=True)
    return compounds_df


def df_to_sdf(compounds_df, file):
    with open(file, 'w') as fi:
        PandasTools.WriteSDF(compounds_df, fi, molColName='ROMol', idName='PID',
                             properties=list(compounds_df.columns))


def df_to_pdb(compounds_df, file):
    with open(file, 'w') as fi:
        pdbwriter = Chem.PDBWriter(fi)
        for _, compound in compounds_df.iterrows():
            mol = Chem.MolFromSmiles(compound.Smiles)
            pdbwriter.write(mol)
        pdbwriter.close()


def df_to_mol(compounds_df, file):
    with open(file, 'w+') as fi:
        for _, compound in compounds_df.iterrows():
            mol = Chem.MolFromSmiles(compound.Smiles)
            molblock = Chem.MolToMolBlock(mol)
            fi.write(molblock)


def update_sdf():
    compounds_df = pd.DataFrame(list(Compound.objects.all().values())).drop('id', axis=1)
    PandasTools.AddMoleculeColumnToFrame(compounds_df, 'Smiles', 'ROMol', includeFingerprints=True)
    df_to_sdf(compounds_df, 'media/all_data.sdf')


import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools


def query_to_df(queryset):
    # Dataframe to write calculations of each compounds
    compounds_df = pd.DataFrame(
        columns=['ID', 'Smiles', 'Molecular_Formula', 'Molecular_Weight', 'H_Bond_Acceptors',
                 'H_Bond_Donors', 'Molar_Refractivity', 'TPSA'])
    for compound in queryset:
        compounds_df = compounds_df.append({
            'ID': compound.PID,
            'Smiles': compound.Smiles,
            'Molecular_Formula': compound.Molecular_Formula,
            'Molecular_Weight': compound.Molecular_Weight,
            'H_Bond_Acceptors': compound.H_Bond_Acceptors,
            'H_Bond_Donors': compound.H_Bond_Donors,
            'Molar_Refractivity': compound.Molar_Refractivity,
            'TPSA': compound.TPSA,
        }, ignore_index=True)
    PandasTools.AddMoleculeColumnToFrame(compounds_df, 'Smiles', 'ROMol', includeFingerprints=True)
    return compounds_df


def df_to_sdf(compounds_df, tmp):
    with open(tmp.name, 'w') as fi:
        PandasTools.WriteSDF(compounds_df, fi, molColName='ROMol', idName='ID',
                             properties=list(compounds_df.columns))


def df_to_pdb(compounds_df, tmp):
    with open(tmp.name, 'w') as fi:
        pdbwriter = Chem.PDBWriter(fi)
        for _, compound in compounds_df.iterrows():
            mol = Chem.MolFromSmiles(compound.Smiles)
            pdbwriter.write(mol)
        pdbwriter.close()


def df_to_mol(compounds_df, tmp):
    with open(tmp.name, 'w+') as fi:
        for _, compound in compounds_df.iterrows():
            mol = Chem.MolFromSmiles(compound.Smiles)
            molblock = Chem.MolToMolBlock(mol)
            fi.write(molblock)


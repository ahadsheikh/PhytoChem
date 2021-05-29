from io import StringIO
import csv

import scrapes

import django
django.setup()

from main.models import Compound

compounds = Compound.objects.all()
c = 1
for compound in compounds[6:8]:
    print("Scrapping data for " + compound.Smiles + "...")
    csvdata = StringIO(scrapes.single_smiles_scrape(compound.Smiles))
    print("Done.")

    csv_reader = csv.reader(csvdata, delimiter=',')

    lc = 0
    for row in csv_reader:
        if lc == 1:
            print("Filling ID: " + compound.PID+ ", Formula:" + compound.Molecular_Formula + "...")

            compound.Molecule_Name = row[0]
            compound.Heavy_Atoms = row[4]
            compound.Arom_Heavy_Atoms = row[5]
            compound.Fraction_Csp3 = row[6]
            compound.Rotatable_bonds = row[7]
            compound.iLOGP = row[12]
            compound.GI_absorption = row[30]
            compound.BBB_permeant = row[31]
            compound.P_gp_substrate = row[32]
            compound.CYP1A2_inhibitor = row[33]
            compound.CYP2C19_inhibitor = row[34]
            compound.CYP2C9_inhibitor = row[35]
            compound.CYP2D6_inhibitor = row[36]
            compound.CYP3A4_inhibitor = row[37]
            compound.LogKp = row[38]
            compound.Lipinski = row[39]
            compound.Ghose = row[40]
            compound.Veber = row[41]
            compound.Egan = row[42]
            compound.Bioavailability_Score = row[44]
            compound.PAINS = row[45]
            compound.Brenk = row[46]
            compound.Leadlikeness = row[47]
            compound.Synthetic_accessibility = row[48]

            compound.save()
            print(f"Done {c} compound\n")
            c += 1

        lc += 1

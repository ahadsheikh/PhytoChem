import csv
from main.models import Compound
import django


def read_csv_to_db(path):
    with open(path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count != 0:
                # Get values from list
                PID = row[1]
                Smiles = row[2]
                Molecular_Formula = row[3]
                Molecular_Weight = row[4]
                H_Bond_Acceptors = row[5]
                H_Bond_Donors = row[6]
                Molar_Refractivity = row[7]
                TPSA = row[8]
                ROMol = row[9]

                # save the values in database
                compound = Compound()
                compound.PID = PID
                compound.Smiles = Smiles
                compound.Molecular_Formula = Molecular_Formula
                compound.Molecular_Weight = Molecular_Weight
                compound.H_Bond_Acceptors = H_Bond_Acceptors
                compound.H_Bond_Donors = H_Bond_Donors
                compound.Molar_Refractivity = Molar_Refractivity
                compound.TPSA = TPSA
                compound.ROMol = ROMol
                compound.save()

                if line_count % 100 == 0:
                    print(line_count, 'Data Processed..')

                # print(f'{id} {PID} {Smiles} {Molecular_Formula} {Molecular_Weight} {H_Bond_Acceptors} {H_Bond_Donors} {Molar_Refractivity} {TPSA}')
                line_count += 1

            else:
                line_count += 1

        print('Processed Successfully.')





import csv
from main.models import Compound


def read_csv_todb():
    with open('utills/data/out.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count != 0 and line_count < 10:
                id = row[0]
                PID = row[1]
                Smiles = row[2]
                Molecular_Formula = row[3]
                Molecular_Weight = row[4]
                H_Bond_Acceptors = row[5]
                H_Bond_Donors = row[6]
                Molar_Refractivity = row[7]
                TPSA = row[8]
                ROMol = row[9]

                print(Smiles)
                line_count += 1
            else:
                line_count += 1

        print(f'Processed {line_count} lines.')





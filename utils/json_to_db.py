import json
from main.models import Compound, Plant


def json_to_db(path):
    with open(path) as json_file:
        data = json.load(json_file)
        data_len = Compound.objects.all().count()
        i = 0
        for k, v in data.items():
            # get data from the dictionary
            Smiles = k
            Molecular_Formula = v['molecular_formula']
            Molecular_Weight = v['molecular_weight']
            H_Bond_Acceptors = v['hba']
            H_Bond_Donors = v['hbd']
            Molar_Refractivity = v['molar_refractivity']
            TPSA = v['tpsa']
            logP = v['logp']
            ROMol = v['romol']
            plants = v['plant']

            # check if this compound already exists
            cur_compound = Compound.objects.filter(Smiles=Smiles)
            if not cur_compound:
                # save the values in database
                compound = Compound()
                compound.PID = 'Phytochem_{}'.format(str(data_len + i).zfill(6))
                compound.Smiles = Smiles
                compound.Molecular_Formula = Molecular_Formula
                compound.Molecular_Weight = Molecular_Weight
                compound.H_Bond_Acceptors = H_Bond_Acceptors
                compound.H_Bond_Donors = H_Bond_Donors
                compound.Molar_Refractivity = Molar_Refractivity
                compound.TPSA = TPSA
                compound.logP = logP
                compound.ROMol = ROMol
                compound.save()
                cur_compound = [compound]
                i += 1
            cur_compound = cur_compound[0]
            # create plants_object list for new plants of this compound
            plants_object = []
            for plant in plants:
                try:
                    # if this plant already exists, just append to the plants_object list
                    current_plant = Plant.objects.get(name=plant)
                    plants_object.append(current_plant)
                except Plant.DoesNotExist:
                    current_plant = Plant.objects.create(name=plant)
                    plants_object.append(current_plant)
            # add all the plants_object
            cur_compound.plants.add(*plants_object)

            if i % 100 == 0:
                print('{} data processed'.format(i))
from django.db import models


class Compound(models.Model):
    PID = models.CharField(max_length=60, unique=True)
    Smiles = models.CharField(max_length=250, db_index=True, unique=True)
    Molecular_Formula = models.CharField(max_length=50, db_index=True)
    Molecular_Weight = models.FloatField()
    H_Bond_Acceptors = models.IntegerField()
    H_Bond_Donors = models.IntegerField()
    Molar_Refractivity = models.FloatField()
    TPSA = models.FloatField()
    ROMol = models.TextField(max_length=150000)

    def __str__(self):
        return self.Molecular_Formula


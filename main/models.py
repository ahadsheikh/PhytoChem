from django.db import models


class Compound(models.Model):
    PID = models.CharField(max_length=60)
    Smiles = models.CharField(max_length=250)
    Molecular_Formula = models.CharField(max_length=50)
    Molecular_Weight = models.FloatField()
    H_Bond_Acceptors = models.IntegerField()
    H_Bond_Donors = models.IntegerField()
    Molar_Refractivity = models.FloatField()
    TPSA = models.FloatField()
    ROMol = models.TextField(max_length=150000)

    def __str__(self):
        return self.Molecular_Formula

from django.db import models


class Plant(models.Model):
    name = models.CharField(max_length=50, unique=True)

    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    def __str__(self):
        return self.name


class Compound(models.Model):
    PID = models.CharField(max_length=60, unique=True)
    Smiles = models.CharField(max_length=5000, db_index=True, unique=True)
    Molecular_Formula = models.CharField(max_length=50, db_index=True)
    Molecular_Weight = models.FloatField()
    H_Bond_Acceptors = models.IntegerField()
    H_Bond_Donors = models.IntegerField()
    Molar_Refractivity = models.FloatField()
    TPSA = models.FloatField()
    ROMol = models.TextField(max_length=150000)
    logP = models.FloatField()
    plants = models.ManyToManyField(Plant, blank=True)

    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    def __str__(self):
        return self.Molecular_Formula


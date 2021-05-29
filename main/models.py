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
    ROMol = models.TextField(max_length=150000)
    Molecule_Name = models.CharField(max_length=200, blank=True)

    # Plant Relation
    plants = models.ManyToManyField(Plant, blank=True)

    # Physicochemical Properties
    Molecular_Formula = models.CharField(max_length=50, db_index=True)
    Molecular_Weight = models.FloatField(blank=True)
    Heavy_Atoms = models.IntegerField(blank=True)
    Arom_Heavy_Atoms = models.IntegerField(blank=True)
    Fraction_Csp3 = models.FloatField(blank=True)
    Rotatable_bonds = models.IntegerField(blank=True)
    H_Bond_Acceptors = models.IntegerField(blank=True)
    H_Bond_Donors = models.IntegerField(blank=True)
    Molar_Refractivity = models.FloatField(blank=True)
    TPSA = models.FloatField(blank=True)

    # Lipophilicity
    iLOGP = models.FloatField(blank=True)

    # Pharmacokinetics
    GI_absorption = models.CharField(max_length=10, blank=True)
    BBB_permeant = models.CharField(max_length=10, blank=True)
    P_gp_substrate = models.CharField(max_length=10, blank=True)
    CYP1A2_inhibitor = models.CharField(max_length=10, blank=True)
    CYP2C19_inhibitor = models.CharField(max_length=10, blank=True)
    CYP2C9_inhibitor = models.CharField(max_length=10, blank=True)
    CYP2D6_inhibitor = models.CharField(max_length=10, blank=True)
    CYP3A4_inhibitor = models.CharField(max_length=10, blank=True)
    LogKp = models.FloatField(blank=True)

    # Druglikeness
    Lipinski = models.IntegerField(blank=True)
    Ghose = models.IntegerField(blank=True)
    Veber = models.IntegerField(blank=True)
    Egan = models.IntegerField(blank=True)
    Muegge = models.IntegerField(blank=True)
    Bioavailability_Score = models.FloatField(blank=True)

    # Medicinal Chemistry
    PAINS = models.CharField(max_length=10, blank=True)
    Brenk = models.CharField(max_length=10, blank=True)
    Leadlikeness = models.IntegerField(blank=True)
    Synthetic_accessibility = models.FloatField(blank=True)

    # Meta
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    def __str__(self):
        return self.Molecular_Formula

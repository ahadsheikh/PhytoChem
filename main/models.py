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
    Molecule_Name = models.CharField(max_length=200, null=True, blank=True)

    # Plant Relation
    plants = models.ManyToManyField(Plant, blank=True)

    # Physicochemical Properties
    Molecular_Formula = models.CharField(max_length=50, db_index=True, null=True, blank=True)
    Molecular_Weight = models.FloatField(null=True, blank=True)
    Heavy_Atoms = models.IntegerField(null=True, blank=True)
    Arom_Heavy_Atoms = models.IntegerField(null=True, blank=True)
    Fraction_Csp3 = models.FloatField(null=True, blank=True)
    Rotatable_bonds = models.IntegerField(null=True, blank=True)
    H_Bond_Acceptors = models.IntegerField(null=True, blank=True)
    H_Bond_Donors = models.IntegerField(null=True, blank=True)
    Molar_Refractivity = models.FloatField(null=True, blank=True)
    TPSA = models.FloatField(null=True, blank=True)

    # Lipophilicity
    iLOGP = models.FloatField(null=True, blank=True)

    # Pharmacokinetics
    GI_absorption = models.CharField(max_length=10, null=True, blank=True)
    BBB_permeant = models.CharField(max_length=10, null=True, blank=True)
    P_gp_substrate = models.CharField(max_length=10, null=True, blank=True)
    CYP1A2_inhibitor = models.CharField(max_length=10, null=True, blank=True)
    CYP2C19_inhibitor = models.CharField(max_length=10, null=True, blank=True)
    CYP2C9_inhibitor = models.CharField(max_length=10, null=True, blank=True)
    CYP2D6_inhibitor = models.CharField(max_length=10, null=True, blank=True)
    CYP3A4_inhibitor = models.CharField(max_length=10, null=True, blank=True)
    LogKp = models.FloatField(null=True, blank=True)

    # Druglikeness
    Lipinski = models.IntegerField(null=True, blank=True)
    Ghose = models.IntegerField(null=True, blank=True)
    Veber = models.IntegerField(null=True, blank=True)
    Egan = models.IntegerField(null=True, blank=True)
    Muegge = models.IntegerField(null=True, blank=True)
    Bioavailability_Score = models.FloatField(null=True, blank=True)

    # Medicinal Chemistry
    PAINS = models.CharField(max_length=10, null=True, blank=True)
    Brenk = models.CharField(max_length=10, null=True, blank=True)
    Leadlikeness = models.IntegerField(null=True, blank=True)
    Synthetic_accessibility = models.FloatField(null=True, blank=True)

    # Meta
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    def __str__(self):
        return self.Molecular_Formula

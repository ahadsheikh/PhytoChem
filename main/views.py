from django.shortcuts import render
from django.db.models import Q
from main.models import Compound


def index(request):
    context = {
        'title': "Home"
    }
    return render(request, 'main/index.html', context=context)


def results(request):
    search = request.GET['search']
    if search.isnumeric():
        search = 'Phytochem_' + search.zfill(5)
    compounds = Compound.objects.filter(Q(PID=search) | Q(Smiles=search) | Q(Molecular_Formula=search))
    context = {
        'title': 'Search results',
        'compounds': compounds
    }
    return render(request, 'main/results.html', context=context)

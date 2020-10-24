from django.shortcuts import render
from main.models import Compound


def index(request):
    context = {
        'title': "Home",
    }
    if request.method == 'POST':
        option = request.POST['option'][0]
        if option == '1':
            pid = request.POST.get('search', -1)
            try:
                compounds = Compound.objects.filter(PID=pid)
                context['compounds'] = compounds
            except Compound.DoesNotExist:
                context['compounds'] = None
        elif option == '2':
            smiles = request.POST.get('search', -1)
            try:
                compounds = Compound.objects.filter(Smiles=smiles)
                context['compounds'] = compounds
            except Compound.DoesNotExist:
                context['compounds'] = None
        else:
            formula = request.POST.get('search', -1)
            try:
                compounds = Compound.objects.filter(Molecular_Formula=formula)
                context['compounds'] = compounds
            except Compound.DoesNotExist:
                context['compounds'] = None

    return render(request, 'main/index.html', context=context)

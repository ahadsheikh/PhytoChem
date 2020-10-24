from django.shortcuts import render
from main.models import Compound

def index(request):
    context = {
        'title': "Home",
    }
    if request.method == 'POST':
        data = request.POST
        if data['option'][0] == '1':
            pid = request.POST.get('search', -1)
            try:
                compounds = Compound.objects.filter(PID=pid)
                context['compounds'] = compounds
            except Compound.DoesNotExist:
                context['compounds'] = None
        elif data['option'][0] == '2':
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
                print(compounds)
                context['compounds'] = compounds
            except Compound.DoesNotExist:
                context['compounds'] = None

        # print(compounds.PID)

    return render(request, 'main/index.html', context=context)

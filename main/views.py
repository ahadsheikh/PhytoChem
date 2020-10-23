from django.shortcuts import render
from main.models import Compound

def index(request):
    context = {
        'title': "Home",
    }
    if request.method == 'POST':
        pid = request.POST.get('search')
        compounds = Compound.objects.get(PID=pid)
        print(compounds.PID)
        context['compounds'] = compounds
    return render(request, 'main/index.html', context=context)

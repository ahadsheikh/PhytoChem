from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
from django.shortcuts import render, redirect
from django.contrib import messages
from dashboard.forms import UploadFileForm


@login_required(redirect_field_name='next')
def dash_index(request):
    context = {
        'title': 'Dashboard'
    }
    return render(request, 'dashboard/dash.html', context=context)


@login_required(redirect_field_name='next')
def upload(request):
    if request.method == 'POST':
        print(request.POST, request.FILES)
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            handle_file(request.FILES['file'])
            messages.success(request, "Upload File Successfully")
            return redirect('dashboard')

    messages.success(request, "Upload File Not Successfully")
    return HttpResponse(
        "<h2>You are not permitted.</h2>"
    )


# Not path view
def handle_file(f):
    with open('media/upload/' + f.__str__(), 'wb') as des:
        for chunk in f.chunks():
            des.write(chunk)

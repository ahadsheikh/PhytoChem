from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
from django.shortcuts import render, redirect
from django.contrib import messages

import os
from dashboard.forms import UploadFileForm
from utils.QueryHandler import handle_new_sdf


@login_required(redirect_field_name='next')
def dash_index(request):
    context = {
        'title': 'Dashboard'
    }
    return render(request, 'dashboard/dash.html', context=context)


@login_required(redirect_field_name='next')
def upload(request):
    if request.method == 'POST':
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
    filename = f.__str__()

    if not os.path.isdir('media/upload'):
        os.mkdir('media')
        os.mkdir('media/upload')
    with open('media/upload/' + filename, 'wb') as des:
        for chunk in f.chunks():
            des.write(chunk)
    try:
        sdf_file = 'media/upload/' + filename
        handle_new_sdf(sdf_file)
    finally:
        os.remove('media/upload/' + filename)

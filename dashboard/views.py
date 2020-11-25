from django.contrib.auth.decorators import login_required
from django.http import HttpResponse, HttpResponseNotFound
from django.shortcuts import render, redirect, get_object_or_404
from django.contrib import messages

import os

from django.utils.decorators import decorator_from_middleware

from core.middlewares import AdminLoginMiddleware
from dashboard.forms import UploadFileForm
from main.views import prepare_download
from submit_data.models import Contribution
from core.utils.QueryHandler import handle_new_sdf


@login_required(redirect_field_name='next')
@decorator_from_middleware(AdminLoginMiddleware)
def dash_index(request):
    contributors = Contribution.objects.all()
    context = {
        'title': 'Dashboard',
        'contributors': contributors
    }
    return render(request, 'dashboard/dash.html', context=context)


@login_required(redirect_field_name='next')
def upload(request):
    if request.method == 'POST' and request.user.is_superuser:
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            handle_file(request.FILES['file'], request.POST['plant'])
            messages.success(request, "Upload File Successfully")
            return redirect('dashboard')

    messages.success(request, "File Upload Failed")
    return HttpResponse(
        "<h2>You are not permitted to view this page.</h2>"
    )


@login_required()
def show_submitted_files(request, cid):
    contribution = get_object_or_404(Contribution, pk=cid)
    new_df = handle_new_sdf(contribution.file.path, change_db=False)
    context = {
        'new_data': new_df.values.tolist(),
        'contributor': contribution.user.first_name + ' ' + contribution.user.last_name,
        'pub_link': contribution.pub_link,
        'data_desc': contribution.data_description,
        'mendeley_data': contribution.mendeley_data_link,
        'plant': contribution.plant_name
    }
    return render(request, 'dashboard/show_data.html', context=context)


@login_required()
def download_new_file(request, cid):
    contribution = get_object_or_404(Contribution, pk=cid)
    return prepare_download(contribution.file.path)


@login_required()
@decorator_from_middleware(AdminLoginMiddleware)
def reject_contribution(request, cid):
    contribution = get_object_or_404(Contribution, pk=cid)
    contribution.status = 2
    contribution.save()
    path = contribution.file.path
    try:
        os.remove(path)
    except FileNotFoundError:
        pass

    return redirect('dashboard')


# Not path view
def handle_file(f, plant=None):
    filename = f.__str__()

    if not os.path.isdir('media/upload'):
        os.mkdir('media')
        os.mkdir('media/upload')
    with open('media/upload/' + filename, 'wb') as des:
        for chunk in f.chunks():
            des.write(chunk)
    try:
        sdf_file = 'media/upload/' + filename
        handle_new_sdf(sdf_file, plant)
    finally:
        os.remove('media/upload/' + filename)

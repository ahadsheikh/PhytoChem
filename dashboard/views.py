from django.http import FileResponse
from django.shortcuts import render, redirect, get_object_or_404
from django.contrib import messages
from django.views import View
from django.views.generic import ListView, DetailView
import os
from core.save_to_db import save_to_db
from dashboard.forms import UploadFileForm
from data_submission.models import Contribution


class SubmissionView(ListView):
    model = Contribution
    template_name = 'dashboard/dash.html'

    def get_queryset(self):
        contributors = Contribution.objects.filter(status=0).order_by('-created_at')
        return contributors


class UploadView(View):
    def post(self, request):
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            handle_file(request.FILES['file'], request.POST['plant'])
            messages.success(request, "Upload File Successfully")
            return redirect('dash:dashboard')
        else:
            messages.error(request, "Data Validation Failed. See below for details.")
            return render(request, 'dashboard/dash.html', {'form': form})

    def get(self, request):
        form = UploadFileForm()
        return render(request, 'dashboard/dash.html', {'form': form})


class SubmittedFileDetailView(DetailView):
    model = Contribution
    template_name = 'dashboard/show_data.html'
    context_object_name = 'contribution'


class SubmittedFileDownloadView(View):
    def get(self, request, cid):
        contribution = Contribution.objects.get(id=cid)
        response = FileResponse(open(contribution.file.path, 'rb'), as_attachment=True, content_type='text/plain')
        return response


class RejectSubmissionView(View):
    def get(self, request, cid):
        contribution = get_object_or_404(Contribution, pk=cid)
        contribution.status = 2
        contribution.save()
        path = contribution.file.path
        try:
            os.remove(path)
        except FileNotFoundError:
            pass
        return redirect('dash:dashboard')


class AcceptSubmissionView(View):
    def get(self, request, cid):
        contribution = get_object_or_404(Contribution, pk=cid)
        contribution.status = 1
        contribution.save()
        path = contribution.file.path
        plant = contribution.plant_name
        save_to_db(plant, path, 0)
        try:
            os.remove(path)
        except FileNotFoundError:
            pass
        return redirect('dash:dashboard')


def handle_file(path, plant):
    filename = path.__str__()
    os.makedirs('media/upload', exist_ok=True)
    with open('media/upload/' + filename, 'wb') as des:
        for chunk in path.chunks():
            des.write(chunk)
    try:
        save_to_db(plant, 'media/upload/' + filename, 0)
    finally:
        os.remove('media/upload/' + filename)

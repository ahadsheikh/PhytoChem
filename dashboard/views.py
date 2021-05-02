from django.http import FileResponse
from django.shortcuts import render, redirect, get_object_or_404
from django.contrib import messages

import os

from django.views import View
from django.views.generic import ListView, DetailView

from dashboard.forms import UploadFileForm
from data_submission.models import Contribution
from core.utils.QueryHandler import handle_new_sdf


class SubmissionView(ListView):
    model = Contribution
    template_name = 'dashboard/dash.html'

    def get_queryset(self):
        contributors = Contribution.objects.all().order_by('-created_at')
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

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        new_df = handle_new_sdf(self.get_object().file.path, change_db=False)
        context['new_data'] = new_df.values.tolist()
        return context


class SubmittedFileDownloadView(View):
    def get(self, request, cid):
        contribution = Contribution.objects.get(id=cid)
        response = FileResponse(open(contribution.file.path, 'rb'), as_attachment=True, content_type='text/plain')
        return response


class RejectSubmissionView(View):
    def post(self, request, cid):
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
    def post(self, request, cid):
        contribution = get_object_or_404(Contribution, pk=cid)
        contribution.status = 1
        contribution.save()
        path = contribution.file.path
        handle_new_sdf(path, plant=contribution.plant_name)
        try:
            os.remove(path)
        except FileNotFoundError:
            pass
        return redirect('dash:dashboard')


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

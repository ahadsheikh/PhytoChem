import os
from django.contrib.auth.mixins import LoginRequiredMixin
from django.shortcuts import render
from django.contrib import messages
from django.shortcuts import redirect
from django.views import View
from data_submission.forms import ContributionForm


class UploadView(LoginRequiredMixin, View):
    def post(self, request):
        form = ContributionForm(request.POST, request.FILES)
        if form.is_valid():
            upload_path = 'media/submittedFiles'
            if not os.path.isdir(upload_path):
                os.makedirs(upload_path)
            messages.success(request, "Upload File Successfully")
            form.save()
            return redirect('data_submission:index')
        else:
            messages.error(request, "Data Validation Failed. See below for details.")
            return render(request, 'data_submission/index.html', {'form': form})

    def get(self, request):
        form = ContributionForm()
        return render(request, 'data_submission/index.html', {'form': form})

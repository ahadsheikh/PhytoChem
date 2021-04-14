import os

from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
from django.shortcuts import render
from django.contrib import messages
from django.shortcuts import redirect

from submit_data.forms import ContributionForm


@login_required
def index(request):
    return render(request, 'submit_data/index.html')


@login_required
def upload(request):
    if request.method == 'POST':
        print(request)
        form = ContributionForm(request.POST, request.FILES)
        if form.is_valid():
            upload_path = 'media/submittedFiles'
            if not os.path.isdir(upload_path):
                os.makedirs(upload_path)
            messages.success(request, "Upload File Successfully")
            form.save()
            return redirect('contribute:submit_data')

        else:
            messages.error(request, "Data Validation Failed. See below for details.")
            return render(request, 'submit_data/index.html', {'form': form})
    else:
        form = ContributionForm()

    return render(request, 'submit_data/index.html', {'form': form})

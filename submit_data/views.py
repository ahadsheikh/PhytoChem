from django.http import HttpResponse
from django.shortcuts import render
from django.contrib import messages
from django.shortcuts import redirect

from submit_data.forms import ContributionForm


def index(request):
    return render(request, 'submit_data/index.html')


def upload(request):
    if request.method == 'POST':
        form = ContributionForm(request.POST, request.FILES)
        if form.is_valid():
            messages.success(request, "Upload File Successfully")
            form.save()
        else:
            messages.success(request, "File Upload Failed")
        return redirect('submit_data')

    return HttpResponse(
        "<h2>You are not permitted to view this page.</h2>"
    )

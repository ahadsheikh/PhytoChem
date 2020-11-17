from django.shortcuts import render


def index(request):
    return render(request, 'submit_data/index.html')


def upload(request):
    pass
    # if request.method == 'POST':
    #     form = UploadFileForm(request.POST, request.FILES)
    #     if form.is_valid():
    #         handle_file(request.FILES['file'], request.POST['plant'])
    #         messages.success(request, "Upload File Successfully")
    #         return redirect('dashboard')
    #
    # messages.success(request, "File Upload Failed")
    # return HttpResponse(
    #     "<h2>You are not permitted to view this page.</h2>"
    # )

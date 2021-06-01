from django.urls import path
from .views import SubmissionView, UploadView, SubmittedFileDetailView,\
    SubmittedFileDownloadView, RejectSubmissionView, AcceptSubmissionView

app_name = "dash"

urlpatterns = [
    path('', SubmissionView.as_view(), name='dashboard'),
    path('upload/', UploadView.as_view(), name='dash_upload'),
    path('show_data/<int:pk>/', SubmittedFileDetailView.as_view(), name='show_data'),
    path('download/<int:cid>/', SubmittedFileDownloadView.as_view(), name='download_new_file'),
    path('reject/<int:cid>/', RejectSubmissionView.as_view(), name='reject_submission_data'),
    path('accept/<int:cid>/', AcceptSubmissionView.as_view(), name='accept_submission_data')
]

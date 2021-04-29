from django.urls import path
from .views import UploadView

app_name = "data_submission"

urlpatterns = [
    path('', UploadView.as_view(), name='index')
]

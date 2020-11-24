from django.urls import path
from . import views

app_name = "contribute"

urlpatterns = [
    path('', views.index, name='submit_data'),
    path('upload/', views.upload, name='user_upload')
]

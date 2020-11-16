from django.urls import path
from . import views

urlpatterns = [
    path('', views.index, name='submit_data'),
]

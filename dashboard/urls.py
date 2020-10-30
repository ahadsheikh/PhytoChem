from django.urls import path
from . import views

urlpatterns = [
    path('', views.dash_index, name='dashboard'),
    path('upload/', views.upload, name='upload')
]

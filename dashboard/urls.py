from django.urls import path
from . import views

app_name = "dash"

urlpatterns = [
    path('', views.dash_index, name='dashboard'),
    path('upload/', views.upload, name='dash_upload'),
    path('show_data/<int:cid>/', views.show_submitted_files, name='show_data'),
    path('download_new/<int:cid>/', views.download_new_file, name='download_new_file'),
    path('reject-contribution/<int:cid>/', views.reject_contribution, name='reject_contribution_data')
]

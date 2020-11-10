from django.urls import path
from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('results/', views.results, name='results'),
    path('download/', views.download_file, name='download_all_results'),
    path('main/plant/<int:id>/', views.plant, name='plant'),
    path('main/compound/<int:id>/', views.compound, name='compound'),
    path('main/download_all/', views.download_all_file, name='all-file-download'),

    # Test

]

from django.urls import path
from django.views.generic import TemplateView

from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('contact/', TemplateView.as_view(template_name='main/contact.html'), name='contact'),
    path('results/', views.results, name='results'),
    path('download/', views.download_file, name='download_all_results'),
    path('download_compound/<str:search>', views.download_file, name='download_compound'),
    path('main/plant/<int:id>/', views.plant, name='plant'),
    path('main/compound/<int:id>/', views.compound, name='compound'),
    path('main/download_all/', views.download_all_file, name='all-file-download'),

    # Test

]

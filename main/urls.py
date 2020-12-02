from django.urls import path
from django.views.generic import TemplateView

from . import views


app_name = "main"

urlpatterns = [
    path('', views.index, name='index'),
    path('contact/', TemplateView.as_view(template_name='main/contact.html'), name='contact'),

    path('about/', TemplateView.as_view(template_name='main/about.html'), name='about'),
    path('about/research-team', TemplateView.as_view(template_name='main/research-team.html'), name='about-research-team'),
    path('about/dev-team', TemplateView.as_view(template_name='main/dev-team.html'), name='about-dev-team'),
    path('about/supervisors', TemplateView.as_view(template_name='main/supervisors.html'), name='about-supervisors'),

    path('results/', views.results, name='results'),
    path('download/', views.download_file, name='download_all_results'),
    path('download_compound/<str:search>', views.download_file, name='download_compound'),
    path('main/plant/<int:id>/', views.plant, name='plant'),
    path('main/compound/<int:id>/', views.compound, name='compound'),
    path('main/download_all/', views.download_all_file, name='all-file-download'),
]

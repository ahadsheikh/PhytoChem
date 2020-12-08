from django.urls import path
from django.views.generic import TemplateView

from . import views


app_name = "main"

urlpatterns = [
    path('', views.index, name='index'),
    path('contact/', TemplateView.as_view(template_name='main/contact.html'), name='contact'),

    path('about/', views.about, name='about'),
    path('members/', views.members, name='members'),
    path('members/research-team', TemplateView.as_view(template_name='main/research-team.html'), name='member-research-team'),
    path('members/dev-team', TemplateView.as_view(template_name='main/dev-team.html'), name='member-dev-team'),
    path('members/supervisors', TemplateView.as_view(template_name='main/supervisors.html'), name='member-supervisors'),

    path('results/', views.results, name='results'),
    path('download/', views.download_file, name='download_all_results'),
    path('download_compound/<str:search>', views.download_file, name='download_compound'),
    path('main/plant/<int:id>/', views.plant, name='plant'),
    path('main/compound/<int:id>/', views.compound, name='compound'),
    path('main/download_all/', views.download_all_file, name='all-file-download'),
]

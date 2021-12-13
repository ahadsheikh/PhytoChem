from django.urls import path
from django.views.generic import TemplateView
from .views import AboutView, ManualView, QueryResultListView, PlantCompoundsListView, CompoundDetailView, FileDownloadView, \
    FullDownloadView

app_name = "main"

urlpatterns = [
    path('', TemplateView.as_view(template_name='main/index.html'), name='index'),
    path('contact/', TemplateView.as_view(template_name='main/contact.html'), name='contact'),
    path('about/', AboutView.as_view(), name='about'),
    path('manual/', ManualView.as_view(), name='manual'),

    path('members/', TemplateView.as_view(template_name='main/members.html'), name='members'),
    path('members/research-team', TemplateView.as_view(template_name='main/research-team.html'),
         name='member-research-team'),
    path('members/dev-team', TemplateView.as_view(template_name='main/dev-team.html'),
         name='member-dev-team'),
    path('members/supervisors', TemplateView.as_view(template_name='main/supervisors.html'),
         name='member-supervisors'),

    path('search', QueryResultListView.as_view(), name='results'),
    path('plant/<int:pk>/', PlantCompoundsListView.as_view(), name='plant'),
    path('phytochem/<int:pk>/', CompoundDetailView.as_view(), name='compound'),

    path('download/<str:search>', FileDownloadView.as_view(), name='download_compound'),
    path('download/', FileDownloadView.as_view(), name='download_all_results'),
    path('download/all/', FullDownloadView.as_view(), name='all-file-download'),
]


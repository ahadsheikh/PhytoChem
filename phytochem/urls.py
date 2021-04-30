from django.contrib import admin
from django.urls import path, include
from django.views.generic import TemplateView

urlpatterns = [
    path('admin/', admin.site.urls),
    path('', include('main.urls', namespace='main')),
    path('dashboard/', include('dashboard.urls', namespace='dash')),
    path('data-submission/', include('data_submission.urls', namespace='data_submission')),
    # Auth System
    path('account/', include('account.urls', namespace='account')),
]

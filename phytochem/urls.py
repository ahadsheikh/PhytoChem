from django.contrib import admin
from django.urls import path, include

urlpatterns = [
    path('admin/', admin.site.urls),
    path('', include('main.urls', namespace='main')),
    path('dashboard/', include('dashboard.urls', namespace='dash')),
    path('data_submission/', include('data_submission.urls', namespace='data_submission')),
    # Auth System
    path('user/', include('account.urls', namespace='account')),
]

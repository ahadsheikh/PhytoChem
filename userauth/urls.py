from django.contrib.auth.views import LoginView
from django.urls import path
from django.views.generic import TemplateView

from . import views

urlpatterns = [
    path('register/', views.register, name='register'),
    path('profile/', views.profile, name='profile'),
    path(
        'login/',
        LoginView.as_view(template_name='userauth/login.html'),
        name='login')
    # Test
]

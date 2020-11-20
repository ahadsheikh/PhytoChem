from django.contrib.auth.views import LoginView, LogoutView
from django.urls import path
from django.views.generic import TemplateView

from . import views

urlpatterns = [
    path('register/', views.register, name='register'),
    path('profile/<str:username>/', views.profile, name='profile'),
    # path('login/', views.login_page, name='login'),
    path(
        'login/',
        LoginView.as_view(template_name='userauth/login.html'),
        name='login'),
    path(
        'logout/',
        LogoutView.as_view(template_name='userauth/logout.html'),
        name='logout'),
    path('profile/edit', views.profile_edit, name='profile_edit'),
    # Test
]

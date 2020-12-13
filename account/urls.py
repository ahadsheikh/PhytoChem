from django.contrib.auth.views import LoginView, LogoutView
from django.urls import path
from django.views.generic import TemplateView

from . import views

app_name = "user"

urlpatterns = [
    path('register/', views.register, name='register'),
    path('verifyemail/', views.verify_email, name='verify_email'),
    path('profile/<int:id>/', views.profile, name='profile'),
    # path('login/', views.login_page, name='login'),
    path(
        'login/',
        LoginView.as_view(template_name='account/login.html'),
        name='login'),
    path('logout/', views.logout, name='logout'),
    path('forgot-password/', views.forgot_password, name='forgot_password'),
    path('profile-edit/', views.profile_edit, name='profile_edit'),
    # Test
]

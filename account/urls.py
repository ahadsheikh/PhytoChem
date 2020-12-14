from django.contrib.auth.views import LoginView, LogoutView
from django.urls import path

from . import views
from .views import ActivateAccount

app_name = "user"

urlpatterns = [
    path('register/', views.register, name='register'),
    path('activate/<uidb64>/<token>/', ActivateAccount.as_view(), name='activate'),
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

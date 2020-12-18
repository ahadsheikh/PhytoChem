from django.contrib.auth.views import LoginView,\
    PasswordResetCompleteView, \
    PasswordResetView, \
    PasswordResetDoneView, \
    PasswordResetConfirmView

from django.urls import path, re_path

from . import views
from .views import ActivateAccount

app_name = "user"

urlpatterns = [
    path('profile/<int:id>/', views.profile, name='profile'),
    path('profile-edit/', views.profile_edit, name='profile_edit'),

    path('register/', views.register, name='register'),
    path('activate/<uidb64>/<token>/', ActivateAccount.as_view(), name='activate'),
    path(
        'login/',
        LoginView.as_view(template_name='account/login.html'),
        name='login'),
    path('logout/', views.logout, name='logout'),

    path('/activate/<uidb64>/<token>/', ActivateAccount.as_view(), name='activate'),
    # Test
]

from django.contrib.auth.views import LoginView, \
    PasswordResetCompleteView, \
    PasswordResetView, \
    PasswordResetDoneView, \
    PasswordResetConfirmView

from django.urls import path, re_path, reverse_lazy

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

    path('password_reset/', PasswordResetView.as_view(template_name='password_change/password_reset_form.html',
                                                      email_template_name='password_change/password_reset_email.html',
                                                      success_url=reverse_lazy('user:password_reset_done')),
         name='password_reset'),
    path('password_reset/<uidb64>/<token>/',
         PasswordResetConfirmView.as_view(template_name='password_change/password_reset_confirm.html',
                                          success_url=reverse_lazy('user:password_reset_complete')),
         name='password_reset_confirm'),
    path('password_reset/done/',
         PasswordResetDoneView.as_view(template_name='password_change/password_reset_done.html'),
         name='password_reset_done'),
    path('password_reset/complete/',
         PasswordResetCompleteView.as_view(template_name='password_change/password_reset_complete.html'),
         name='password_reset_complete'),

    path('activate/<uidb64>/<token>/', ActivateAccount.as_view(), name='activate'),
    # Test
]

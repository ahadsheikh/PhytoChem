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

    re_path(r'^password_reset/$',
            PasswordResetView.as_view(template_name='password_change/password_reset_form.html'),
            name='password_reset'),
    re_path(r'^password_reset/done/$',
            PasswordResetDoneView.as_view(template_name='password_change/password_reset_done.html'),
            name='password_reset_done'),
    re_path(
        r'^reset/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$',
        PasswordResetConfirmView.as_view(template_name='password_change/password_reset_confirm.html'),
        name='password_reset_confirm'),
    re_path(r'^reset/done/$',
            PasswordResetCompleteView.as_view(template_name='password_change/password_reset_complete.html'),
            name='password_reset_complete'),
    # Test
]

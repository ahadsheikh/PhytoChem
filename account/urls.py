from django.contrib.auth.views import LoginView, PasswordResetCompleteView, PasswordResetView, \
    PasswordResetDoneView, PasswordResetConfirmView, LogoutView
from django.urls import path, reverse_lazy
from .views import ActivateAccount, PasswordChange, RegisterView, ProfileView, ProfileEditView

app_name = "account"

urlpatterns = [
    path('<int:pk>/', ProfileView.as_view(), name='profile'),
    path('edit_profile/', ProfileEditView.as_view(), name='profile_edit'),

    path('login/', LoginView.as_view(template_name='account/login.html'), name='login'),
    path('logout/', LogoutView.as_view(next_page='main:index'), name='logout'),
    path('register/', RegisterView.as_view(), name='register'),
    path('activate/<uidb64>/<token>/', ActivateAccount.as_view(), name='activate'),

    path('password_reset/',
         PasswordResetView.as_view(template_name='password_change/password_reset_form.html',
                                   email_template_name='password_change/password_reset_email.html',
                                   success_url=reverse_lazy('account:password_reset_done')),
         name='password_reset'),
    path('password_reset/<uidb64>/<token>/',
         PasswordResetConfirmView.as_view(template_name='password_change/password_reset_confirm.html',
                                          success_url=reverse_lazy('account:password_reset_complete')),
         name='password_reset_confirm'),
    path('password_reset/done/',
         PasswordResetDoneView.as_view(template_name='password_change/password_reset_done.html'),
         name='password_reset_done'),
    path('password_reset/complete/',
         PasswordResetCompleteView.as_view(template_name='password_change/password_reset_complete.html'),
         name='password_reset_complete'),
    path('password_change/', PasswordChange.as_view(), name='password_change'),
]

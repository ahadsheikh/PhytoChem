from django.contrib.auth.mixins import LoginRequiredMixin
from django.contrib.sites.shortcuts import get_current_site
from django.contrib.auth import authenticate, login
from django.shortcuts import render, redirect
from django.contrib import messages
from django.contrib.auth.models import User
from django.template.loader import render_to_string
from django.utils.encoding import force_text, force_bytes
from django.utils.http import urlsafe_base64_decode, urlsafe_base64_encode
from django.views.generic import View, DetailView

from data_submission.models import Contribution
from account.forms import AccountForm, AccountUpdateForm, PasswordChangeForm
from .models import Account
from account.tokens import account_activation_token


class RegisterView(View):
    def post(self, request):
        form = AccountForm(request.POST)
        if form.is_valid():
            email = form.cleaned_data['email']
            password = form.cleaned_data['password1']
            form.save()
            user = authenticate(username=email, password=password)
            current_site = get_current_site(request)

            subject = 'Activate Your Phytochem Database Account'
            message = render_to_string('email_verification/email_verification.html', {
                'account': user,
                'domain': current_site.domain,
                'uid': urlsafe_base64_encode(force_bytes(user.pk)),
                'token': account_activation_token.make_token(user),
            })
            user.is_active = False
            user.save()
            if user.email_user(subject, message) == 1:
                messages.success(request, 'Please check your email and confirm the link to complete registration.')
            else:
                messages.warning(request, 'Failed to confirm email')
            return redirect('account:register')
        else:
            return render(request, 'account/register.html', {'form': form})

    def get(self, request):
        form = AccountForm()
        return render(request, 'account/register.html', {'form': form})


class ActivateAccount(View):
    def get(self, request, uidb64, token, *args, **kwargs):
        invalid_link = '''
        The email confirmation link was invalid, possibly because it has already been used.
        Please request a new password reset.
        '''
        expired_link = '''
        User already confirmed registration! Link expired.  
        '''
        try:
            uid = force_text(urlsafe_base64_decode(uidb64))
            user = Account.objects.get(pk=uid)
        except (TypeError, ValueError, OverflowError, User.DoesNotExist):
            user = None
        if user is not None and account_activation_token.check_token(user, token):
            user.is_active = True
            user.save()
            login(request, user)
            messages.success(request, 'Your account have been confirmed.')
            return redirect('account:profile', user.id)
        else:
            return render(request, 'email_verification/invalid_link.html', {'vallidation_error': invalid_link})


class PasswordChange(LoginRequiredMixin, View):
    def get(self, request):
        form = PasswordChangeForm()
        return render(request, 'password_change/password_reset_form.html', {'form': form})

    def post(self, request):
        form = PasswordChangeForm(request.POST)
        if form.is_valid():
            user = request.user
            current_site = get_current_site(request)

            subject = 'Activate Your Phytochem Database Account'
            message = render_to_string('email_verification/email_verification.html', {
                'account': user,
                'domain': current_site.domain,
                'uid': urlsafe_base64_encode(force_bytes(user.pk)),
                'token': account_activation_token.make_token(user),
            })
            if user.email_user(subject, message) == 1:
                return render(request, 'password_change/password_reset_done.html')
            else:
                messages.warning(request, 'Failed to confirm email')
                return render(request, 'password_change/password_reset_form.html', {'form': form})


class ProfileView(DetailView):
    model = Account
    context_object_name = 'user_d'
    template_name = 'account/profile.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        contributions = Contribution.objects.filter(user=self.get_object()).order_by('-created_at')
        context['contributions'] = contributions
        return context


class ProfileEditView(LoginRequiredMixin, View):
    def post(self, request):
        acc_up_form = AccountUpdateForm(request.POST, instance=request.user)
        if acc_up_form.is_valid():
            acc_up_form.save()
            messages.success(request, "Successfully Saved")
            return redirect('account:profile', request.user.id)
        else:
            return render(request, 'account/profile_edit.html', {'form': acc_up_form})

    def get(self, request):
        acc_up_form = AccountUpdateForm(instance=request.user)
        return render(request, 'account/profile_edit.html', {'form': acc_up_form})

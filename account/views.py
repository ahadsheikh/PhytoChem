from django.contrib.sites.shortcuts import get_current_site
from django.contrib.auth import authenticate
from django.contrib.auth.decorators import login_required
from django.shortcuts import render, redirect, get_object_or_404
from django.contrib import messages, auth
from django.contrib.auth import login
from django.contrib.auth.models import User
from django.template.loader import render_to_string
from django.utils.encoding import force_text, force_bytes
from django.utils.http import urlsafe_base64_decode, urlsafe_base64_encode
from django.views.generic import View

from submit_data.models import Contribution
from account.forms import AccountForm, AccountUpdateForm
from .models import Account
from account.tokens import account_activation_token


def logout(request):
    auth.logout(request)
    return redirect('main:index')


def register(request):
    if request.method == 'POST':
        form = AccountForm(request.POST)
        if form.is_valid():
            email = form.cleaned_data['email']
            password = form.cleaned_data['password1']
            form.save()
            user = authenticate(username=email, password=password)
            current_site = get_current_site(request)

            subject = 'Activate Your Phytochem Database Account'
            message = render_to_string('email_verification/email_verification.html', {
                'user': user,
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
            return redirect('user:register')
        else:
            return render(request, 'account/register.html', {'form': form})
    else:
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
            return redirect('user:profile', user.id)
        else:
            # messages.warning(request, 'The confirmation link was invalid, possibly because it has already been used.')
            return render(request, 'email_verification/invalid_link.html', {'vallidation_error': invalid_link})


@login_required
def profile(request, id):
    user = get_object_or_404(Account, id=id)
    contributions = Contribution.objects.filter(user=user).order_by('-created_at')
    context = {
        'user_d': user,
        'contributions': contributions
    }
    return render(request, 'account/profile.html', context=context)


@login_required
def profile_edit(request):
    if request.method == 'POST':
        accUpForm = AccountUpdateForm(request.POST, instance=request.user)
        if accUpForm.is_valid():
            accUpForm.save()
            messages.success(request, "Successfully Saved")
            return redirect('user:profile', request.user.id)
        else:
            return render(request, 'account/profile_edit.html', {'form': accUpForm})

    accUpForm = AccountUpdateForm(instance=request.user)

    return render(request, 'account/profile_edit.html', {'form': accUpForm})

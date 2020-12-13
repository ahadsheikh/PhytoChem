from django.http import HttpResponse
from django.shortcuts import render

from django.contrib.auth import authenticate, login
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.shortcuts import render, redirect, get_object_or_404
from django.contrib import messages, auth

from submit_data.models import Contribution
from account.forms import AccountForm, AccountUpdateForm, ForgotPasswordForm, AcceptForgotPasswordCodeForm
from .models import Account, ForgotPasswordCode


from random import randint


def logout(request):
    auth.logout(request)
    return redirect('main:index')


def forgot_password(request):
    if request.method == 'POST':
        form = ForgotPasswordForm(request.POST)
        if form.is_valid():
            try:
                user = Account.objects.get(email=form.cleaned_data['email'])
                code = randint(100000, 999999)
                forgetCode = ForgotPasswordCode.objects.filter(user=user)
                if not forgetCode:
                    ForgotPasswordCode.objects.create(user=user, code=code)
                return render(request, 'account/submit_code.html', {'email': user.email})
            except Account.DoesNotExist:
                errors = ("Account doesn't exists", )
                return render(request, 'account/forgot_password.html', {'errors': errors})
    return render(request, 'account/forgot_password.html', {})


def code_accept(request):
    if request.method == 'POST':
        accForgotForm = AcceptForgotPasswordCodeForm(request.POST)
        if accForgotForm.is_valid():
            user = None
            try:
                user = Account.objects.get(email=accForgotForm.cleaned_data['email'])
            except Account.DoesNotExist:
                errors = ('Not Valid Code',)
                return render(request, 'account/forgot_password.html', {'errors': errors})

        else:
            errors = ('Not Valid Code', )
            return render(request, 'account/forgot_password.html', {'errors': errors})

    return HttpResponse('Not permitted Request')


def register(request):
    if request.method == 'POST':
        form = AccountForm(request.POST)
        if form.is_valid():
            email = form.cleaned_data['email']
            password = form.cleaned_data['password1']
            form.save()
            user = authenticate(username=email, password=password)
            login(request, user)
            return redirect('user:verify_email')
        else:
            return render(request, 'account/register.html', {'form': form})
    else:
        form = AccountForm()
    return render(request, 'account/register.html', {'form': form})


def verify_email(request):
    return HttpResponse("Hi")


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




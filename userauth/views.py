from django.contrib.auth import authenticate, login
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.shortcuts import render, redirect, get_object_or_404
from django.contrib import messages

from submit_data.models import Contribution
from userauth.forms import UserForm, UserUpdateForm, ProfileUpdateForm, LoginUsernameForm, LoginEmailForm
from .models import Profile


def register(request):
    if request.method == 'POST':
        form = UserForm(request.POST)
        if form.is_valid():
            username = form.cleaned_data['username']
            password = form.cleaned_data['password1']
            form.save()
            user = authenticate(username=username, password=password)
            login(request, user)
            return redirect('user:profile', username)
    else:
        form = UserForm()
    return render(request, 'userauth/register.html', {'form': form})


@login_required
def profile(request, username):
    user = get_object_or_404(User, username=username)
    contributions = Contribution.objects.filter(user=user)
    context = {
        'user_d': user,
        'contributions': contributions
    }
    return render(request, 'userauth/profile.html', context=context)


@login_required
def profile_edit(request):
    if request.method == 'POST':
        userUpdateForm = UserUpdateForm(request.POST, instance=request.user)
        profileUpdateForm = ProfileUpdateForm(request.POST, instance=request.user.profile)

        if userUpdateForm.is_valid() and profileUpdateForm.is_valid():
            userUpdateForm.save()
            profileUpdateForm.save()
            messages.success(request, "Successfully Edit Done")
            return redirect('user:profile_edit')

    userUpdateForm = UserUpdateForm(instance=request.user)
    profileUpdateForm = ProfileUpdateForm(instance=request.user.profile)

    context = {
        'u_form': userUpdateForm,
        'p_form': profileUpdateForm
    }

    return render(request, 'userauth/profile_edit.html', context=context)



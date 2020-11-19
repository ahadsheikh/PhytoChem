from django.contrib.auth import authenticate, login
from django.contrib.auth.decorators import login_required
from django.shortcuts import render, redirect
from django.contrib import messages

from userauth.forms import UserForm, UserUpdateForm, ProfileUpdateForm


def register(request):
    context = {
        'title': 'Register',
    }
    if request.method == 'POST':
        form = UserForm(request.POST)
        if form.is_valid():
            username = form.cleaned_data['username']
            password = form.cleaned_data['password1']
            form.save()
            user = authenticate(username=username, password=password)
            login(request, user)
            return redirect('profile')
    else:
        form = UserForm()
    return render(request, 'userauth/register.html', {'form': form})


def profile(request):
    return render(request, 'userauth/profile.html', {})


@login_required
def profile_edit(request):
    if request.method == 'POST':
        userUpdateForm = UserUpdateForm(request.POST, instance=request.user)
        profileUpdateForm = ProfileUpdateForm(request.POST, instance=request.user.profile)

        if userUpdateForm.is_valid() and profileUpdateForm.is_valid():
            userUpdateForm.save()
            profileUpdateForm.save()
            messages.success(request, "Successfully Edit Done")
            return redirect('profile_edit')

    userUpdateForm = UserUpdateForm(instance=request.user)
    profileUpdateForm = ProfileUpdateForm(instance=request.user.profile)

    context = {
        'title': 'Profile | Edit',
        'u_form': userUpdateForm,
        'p_form': profileUpdateForm
    }

    return render(request, 'userauth/profile_edit.html', context=context)



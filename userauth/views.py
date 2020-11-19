from django.contrib.auth import authenticate, login
from django.shortcuts import render, redirect

from userauth.forms import UserForm


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

from django.http import HttpResponse
from django.shortcuts import render

from django.contrib.auth import authenticate, login
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.shortcuts import render, redirect, get_object_or_404
from django.contrib import messages

from submit_data.models import Contribution
from account.forms import AccountForm, AccountUpdateForm
from .models import Account


def register(request):
    if request.method == 'POST':
        form = AccountForm(request.POST)
        if form.is_valid():
            email = form.cleaned_data['email']
            password = form.cleaned_data['password1']
            form.save()
            user = authenticate(username=email, password=password)
            login(request, user)
            return redirect('user:profile', user.id)
        else:
            return render(request, 'account/register.html', {'form': form})
    else:
        form = AccountForm()
    return render(request, 'account/register.html', {'form': form})


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
        print(request.POST)
        if accUpForm.is_valid():
            print(accUpForm.cleaned_data)
            accUpForm.save()
            messages.success(request, "Successfully Edit Done")
            return redirect('user:profile_edit')

    accUpForm = AccountUpdateForm(instance=request.user)

    return render(request, 'account/profile_edit.html', {'form': accUpForm})




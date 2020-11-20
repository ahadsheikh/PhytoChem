from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User
from django import forms

from .models import Profile


class UserForm(UserCreationForm):
    username = forms.CharField(max_length=100)
    first_name = forms.CharField(max_length=20)
    last_name = forms.CharField(max_length=20)
    email = forms.EmailField(max_length=200)

    class Meta:
        model = User
        fields = ['username', 'email', 'first_name', 'last_name', 'password1', 'password2']


class UserUpdateForm(forms.ModelForm):
    class Meta:
        model = User
        fields = ['username', 'email', 'first_name', 'last_name']


class ProfileUpdateForm(forms.ModelForm):

    class Meta:
        model = Profile
        fields = ['about']


class LoginUsernameForm(forms.Form):
    uname = forms.CharField(max_length=100)
    password = forms.PasswordInput()


class LoginEmailForm(forms.Form):
    uname = forms.EmailField(max_length=200)
    password = forms.PasswordInput()

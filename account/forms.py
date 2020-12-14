from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User
from django import forms
from django.core.validators import MinValueValidator

from .models import Account


class AccountForm(UserCreationForm):
    first_name = forms.CharField(max_length=20)
    last_name = forms.CharField(max_length=20)

    class Meta:
        model = Account
        fields = ['email', 'first_name', 'last_name', 'password1', 'password2']


class AccountUpdateForm(forms.ModelForm):

    class Meta:
        model = Account
        fields = ['first_name', 'last_name', 'location', 'linkedin', 'researchgate', 'about']

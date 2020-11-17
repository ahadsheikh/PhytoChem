from django import forms

from submit_data.models import Contributor


class ContributorForm(forms.ModelForm):

    class Meta:
        model = Contributor
        fields = ['name', 'email', 'plantName', 'file']
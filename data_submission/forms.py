from django import forms
from data_submission.models import Contribution


class ContributionForm(forms.ModelForm):
    class Meta:
        model = Contribution
        fields = ['user', 'plant_name', 'pub_link', 'data_description', 'mendeley_data_link', 'file']

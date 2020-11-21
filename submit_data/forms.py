from django import forms

from submit_data.models import Contribution


class ContributionForm(forms.ModelForm):

    class Meta:
        model = Contribution
        fields = ['plant_name', 'pub_link', 'data_description', 'mendeley_data_link', 'file']

from django import forms
from django.core.validators import FileExtensionValidator, URLValidator

from core.utils.QueryHandler import validate_sdf
from submit_data.models import Contribution


class ContributionForm(forms.ModelForm):

    class Meta:
        model = Contribution
        fields = ['user', 'plant_name', 'pub_link', 'data_description', 'mendeley_data_link', 'file']

from django import forms
from django.core.validators import FileExtensionValidator

from core.utils.QueryHandler import validate_sdf
from submit_data.models import Contribution


class ContributionForm(forms.ModelForm):

    file = forms.FileField(validators=[FileExtensionValidator(allowed_extensions=['sdf']),
                                        validate_sdf])

    class Meta:
        model = Contribution
        fields = ['user', 'plant_name', 'pub_link', 'data_description', 'mendeley_data_link', 'file']

from django import forms
from django.core.validators import FileExtensionValidator

from core.utils.QueryHandler import validate_sdf


class UploadFileForm(forms.Form):
    file = forms.FileField(validators=[FileExtensionValidator(allowed_extensions=['sdf']),
                                       validate_sdf])

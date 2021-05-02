from django.contrib.auth import get_user_model
from django.core.validators import FileExtensionValidator
from django.db import models
from core.utils.QueryHandler import validate_sdf

user = get_user_model()


class Contribution(models.Model):
    user = models.ForeignKey(user, on_delete=models.CASCADE)
    plant_name = models.CharField(max_length=100)
    pub_link = models.URLField(max_length=200)
    data_description = models.TextField(max_length=500, blank=True)
    mendeley_data_link = models.URLField(max_length=200, blank=True)
    status = models.SmallIntegerField(default=0)  # 0 -> Not reviewed, 1 -> Accepted, 2 -> Rejected
    file = models.FileField(upload_to='submittedFiles',
                            validators=[FileExtensionValidator(allowed_extensions=['sdf']),
                                        validate_sdf])
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.user.email

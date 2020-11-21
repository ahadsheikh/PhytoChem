from django.db import models
from django.contrib.auth.models import User


class Contribution(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    plant_name = models.CharField(max_length=100)
    pub_link = models.URLField(max_length=200)
    data_description = models.TextField(max_length=500, blank=True)
    mendeley_data_link = models.URLField(max_length=200, blank=True)
    file = models.FileField(upload_to='submittedFiles')

    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.name
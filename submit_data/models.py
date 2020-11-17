from django.db import models


class Contributor(models.Model):
    name = models.CharField(max_length=50)
    email = models.EmailField(max_length=50)
    plantName = models.CharField(max_length=100)
    file = models.FileField(upload_to='submittedFiles')

    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.name


class ContributorDatabase(models.Model):
    name = models.CharField(max_length=50)
    email = models.EmailField(max_length=50)

    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.name
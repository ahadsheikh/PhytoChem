# Generated by Django 3.1.2 on 2020-12-02 08:45

from django.conf import settings
import django.core.validators
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='Contribution',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('plant_name', models.CharField(max_length=100)),
                ('pub_link', models.URLField()),
                ('data_description', models.TextField(blank=True, max_length=500)),
                ('mendeley_data_link', models.URLField(blank=True)),
                ('status', models.SmallIntegerField(default=0)),
                ('file', models.FileField(upload_to='submittedFiles', validators=[django.core.validators.FileExtensionValidator(allowed_extensions=['sdf'])])),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('user', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL)),
            ],
        ),
    ]

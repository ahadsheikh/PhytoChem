# Generated by Django 3.1.2 on 2020-10-26 07:42

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('main', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='compound',
            name='Smiles',
            field=models.CharField(db_index=True, max_length=250, unique=True),
        ),
    ]
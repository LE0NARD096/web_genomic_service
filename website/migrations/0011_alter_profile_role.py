# Generated by Django 5.0.1 on 2024-02-09 14:10

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('website', '0010_rename_annotationstatus_annotationstatu'),
    ]

    operations = [
        migrations.AlterField(
            model_name='profile',
            name='role',
            field=models.CharField(choices=[('user', 'user'), ('annotator', 'annotator'), ('validator', 'validator'), ('admin', 'admin')], max_length=50, verbose_name='role'),
        ),
    ]
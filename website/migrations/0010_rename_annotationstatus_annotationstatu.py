# Generated by Django 5.0.1 on 2024-02-09 13:11

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('contenttypes', '0002_remove_content_type_name'),
        ('website', '0009_rename_annotationvalidation_annotationstatus'),
    ]

    operations = [
        migrations.RenameModel(
            old_name='AnnotationStatus',
            new_name='AnnotationStatu',
        ),
    ]

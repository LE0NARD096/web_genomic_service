# Generated by Django 4.2.7 on 2024-01-28 13:36

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('website', '0001_initial'),
    ]

    operations = [
        migrations.RenameField(
            model_name='annotationprotein',
            old_name='annotations',
            new_name='geneprotein',
        ),
    ]

# Generated by Django 5.0.1 on 2024-02-09 11:34

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('website', '0006_remove_annotationgenome_is_validated_and_more'),
    ]

    operations = [
        migrations.AlterField(
            model_name='genome',
            name='chromosome',
            field=models.CharField(max_length=255, unique=True),
        ),
    ]

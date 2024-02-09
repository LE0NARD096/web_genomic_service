# Generated by Django 5.0.1 on 2024-02-09 13:08

import django.db.models.deletion
from django.conf import settings
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('contenttypes', '0002_remove_content_type_name'),
        ('website', '0007_alter_genome_chromosome'),
    ]

    operations = [
        migrations.AlterField(
            model_name='annotationprotein',
            name='geneprotein',
            field=models.OneToOneField(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='website.geneprotein'),
        ),
        migrations.AlterField(
            model_name='profile',
            name='role',
            field=models.CharField(choices=[('user', 'user'), ('annotator', 'annotator'), ('validator', 'validator')], max_length=50, verbose_name='role'),
        ),
        migrations.CreateModel(
            name='AnnotationValidation',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('object_id', models.PositiveIntegerField()),
                ('comment', models.TextField(blank=True)),
                ('validation_time', models.DateTimeField(auto_now=True)),
                ('content_type', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='contenttypes.contenttype')),
                ('validator', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL)),
            ],
        ),
    ]
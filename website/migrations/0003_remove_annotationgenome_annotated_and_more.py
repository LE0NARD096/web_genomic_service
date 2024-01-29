# Generated by Django 4.2.7 on 2024-01-28 14:11

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('website', '0002_rename_annotations_annotationprotein_geneprotein'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='annotationgenome',
            name='annotated',
        ),
        migrations.RemoveField(
            model_name='annotationprotein',
            name='annotated',
        ),
        migrations.AddField(
            model_name='geneprotein',
            name='annotated',
            field=models.BooleanField(default=False, verbose_name='annotated'),
        ),
        migrations.AddField(
            model_name='genome',
            name='annotated',
            field=models.BooleanField(default=False, verbose_name='annotated'),
        ),
        migrations.AlterField(
            model_name='genome',
            name='sequence',
            field=models.TextField(blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='genome',
            name='upload_time',
            field=models.DateTimeField(null=True),
        ),
    ]

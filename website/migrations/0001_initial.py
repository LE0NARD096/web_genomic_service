from django.db import migrations, models


class Migration(migrations.Migration):
    initial = True

    dependencies = []

    operations = [
        migrations.CreateModel(
            name="Genome",
            fields=[
                (
                    "id",
                    models.BigAutoField(
                        auto_created=True,
                        primary_key=True,
                        serialize=False,
                        verbose_name="ID",
                    ),
                ),
                ("sequence", models.TextField()),
                ("species", models.CharField(max_length=255)),
                ("description", models.CharField(max_length=255)),
                ("type", models.CharField(max_length=255)),
            ],
        ),
        migrations.CreateModel(
            name="User",
            fields=[
                (
                    "id",
                    models.BigAutoField(
                        auto_created=True,
                        primary_key=True,
                        serialize=False,
                        verbose_name="ID",
                    ),
                ),
                ("email", models.EmailField(max_length=200)),
                ("password", models.CharField(max_length=200)),
                ("firstName", models.CharField(max_length=200)),
                ("lastName", models.CharField(max_length=200)),
                ("phoneNumber", models.IntegerField()),
            ],
        ),
    ]

from django.db import migrations, models


class Migration(migrations.Migration):
    initial = True

    dependencies = []

    operations = [
        migrations.CreateModel(
            name="Genome",
            fields=[
                (
                    "id",
                    models.BigAutoField(
                        auto_created=True,
                        primary_key=True,
                        serialize=False,
                        verbose_name="ID",
                    ),
                ),
                ("sequence", models.TextField()),
                ("species", models.CharField(max_length=255)),
                ("description", models.CharField(max_length=255)),
                ("type", models.CharField(max_length=255)),
            ],
        ),
        migrations.CreateModel(
            name="User",
            fields=[
                (
                    "id",
                    models.BigAutoField(
                        auto_created=True,
                        primary_key=True,
                        serialize=False,
                        verbose_name="ID",
                    ),
                ),
                ("email", models.EmailField(max_length=200)),
                ("password", models.CharField(max_length=200)),
                ("firstName", models.CharField(max_length=200)),
                ("lastName", models.CharField(max_length=200)),
                ("phoneNumber", models.IntegerField()),
            ],
        ),
    ]

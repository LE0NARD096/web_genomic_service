# Generated by Django 4.2.7 on 2024-01-28 13:33

from django.conf import settings
import django.contrib.auth.models
import django.contrib.auth.validators
from django.db import migrations, models
import django.db.models.deletion
import django.utils.timezone


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        ('auth', '0012_alter_user_first_name_max_length'),
    ]

    operations = [
        migrations.CreateModel(
            name='Profile',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('password', models.CharField(max_length=128, verbose_name='password')),
                ('last_login', models.DateTimeField(blank=True, null=True, verbose_name='last login')),
                ('is_superuser', models.BooleanField(default=False, help_text='Designates that this user has all permissions without explicitly assigning them.', verbose_name='superuser status')),
                ('username', models.CharField(error_messages={'unique': 'A user with that username already exists.'}, help_text='Required. 150 characters or fewer. Letters, digits and @/./+/-/_ only.', max_length=150, unique=True, validators=[django.contrib.auth.validators.UnicodeUsernameValidator()], verbose_name='username')),
                ('is_staff', models.BooleanField(default=False, help_text='Designates whether the user can log into this admin site.', verbose_name='staff status')),
                ('is_active', models.BooleanField(default=True, help_text='Designates whether this user should be treated as active. Unselect this instead of deleting accounts.', verbose_name='active')),
                ('date_joined', models.DateTimeField(default=django.utils.timezone.now, verbose_name='date joined')),
                ('email', models.EmailField(max_length=200, unique=True, verbose_name='email')),
                ('phoneNumber', models.CharField(max_length=15, verbose_name='phone_number')),
                ('first_name', models.CharField(max_length=200, verbose_name='first_name')),
                ('last_name', models.CharField(max_length=200, verbose_name='last_name')),
                ('role', models.CharField(max_length=50, verbose_name='role')),
                ('is_approved', models.BooleanField(default=False, verbose_name='approved')),
                ('groups', models.ManyToManyField(blank=True, help_text='The groups this user belongs to. A user will get all permissions granted to each of their groups.', related_name='user_set', related_query_name='user', to='auth.group', verbose_name='groups')),
                ('user_permissions', models.ManyToManyField(blank=True, help_text='Specific permissions for this user.', related_name='user_set', related_query_name='user', to='auth.permission', verbose_name='user permissions')),
            ],
            options={
                'verbose_name': 'user',
                'verbose_name_plural': 'users',
                'abstract': False,
            },
            managers=[
                ('objects', django.contrib.auth.models.UserManager()),
            ],
        ),
        migrations.CreateModel(
            name='Genome',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('sequence', models.TextField()),
                ('chromosome', models.CharField(max_length=255)),
                ('start', models.IntegerField(null=True)),
                ('end', models.IntegerField(null=True)),
                ('upload_time', models.DateTimeField()),
            ],
        ),
        migrations.CreateModel(
            name='GeneProtein',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('accession_number', models.CharField(max_length=255)),
                ('sequence', models.TextField()),
                ('type', models.CharField(max_length=255)),
                ('start', models.IntegerField()),
                ('end', models.IntegerField()),
                ('upload_time', models.DateTimeField()),
                ('genome', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='website.genome')),
            ],
        ),
        migrations.CreateModel(
            name='AnnotationProtein',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('gene', models.CharField(default='unknown', max_length=255)),
                ('transcript', models.CharField(default='unknown', max_length=255)),
                ('gene_biotype', models.CharField(default='unknown', max_length=255)),
                ('transcript_biotype', models.CharField(default='unknown', max_length=255)),
                ('gene_symbol', models.CharField(default='unknown', max_length=255)),
                ('description', models.CharField(default='empty', max_length=1000)),
                ('annotation_time', models.DateTimeField(null=True)),
                ('annotated', models.BooleanField(default=False, verbose_name='annotated')),
                ('annotations', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='website.geneprotein')),
                ('annotator', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL)),
            ],
        ),
        migrations.CreateModel(
            name='AnnotationGenome',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('species', models.CharField(default='unknown', max_length=255)),
                ('annotation_time', models.DateTimeField(null=True)),
                ('annotated', models.BooleanField(default=False, verbose_name='annotated')),
                ('annotator', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL)),
                ('genome', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, to='website.genome')),
            ],
        ),
    ]

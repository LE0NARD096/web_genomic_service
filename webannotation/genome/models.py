from django.db import models
from django.contrib.auth.models import User
from django.contrib.auth.models import AbstractUser
from django.conf import settings
from django.contrib.auth.models import AbstractUser, Group, Permission

class CustomUser(AbstractUser):
    last_connection = models.DateTimeField(null=True, blank=True)
    is_approved = models.BooleanField(default=False)
    role = models.ForeignKey('Role', on_delete=models.SET_NULL, null=True, blank=True)

    groups = models.ManyToManyField(
        Group,
        verbose_name=('groups'),
        blank=True,
        help_text=(
            'The groups this user belongs to. A user will get all permissions '
            'granted to each of their groups.'
        ),
        related_name="customuser_set",
        related_query_name="customuser",
    )
    user_permissions = models.ManyToManyField(
        Permission,
        verbose_name=('user permissions'),
        blank=True,
        help_text=('Specific permissions for this user.'),
        related_name="customuser_set",
        related_query_name="customuser",
    )

    class Meta:
        db_table = 'custom_user'



class Role(models.Model):
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

class Genome(models.Model):
    name = models.CharField(max_length=200)

    def __str__(self):
        return self.name

class File(models.Model):
    filePath = models.CharField(max_length=200)

    def __str__(self):
        return self.filePath

class Sequence(models.Model):
    isUnique = models.BooleanField(default=False)
    genome = models.ForeignKey(Genome, on_delete=models.CASCADE)
    cdsFile = models.ForeignKey(File, related_name='cds_file', on_delete=models.CASCADE)
    pepFile = models.ForeignKey(File, related_name='pep_file', on_delete=models.CASCADE)

    def __str__(self):
        return f"Sequence {self.id}"



class Annotation(models.Model):
    text = models.TextField()
    user = models.ForeignKey(settings.AUTH_USER_MODEL, on_delete=models.CASCADE)
    sequence = models.ForeignKey(Sequence, on_delete=models.CASCADE)
    isValidated = models.BooleanField(default=False)
    createdAt = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"Annotation {self.id} by {self.user.username}"


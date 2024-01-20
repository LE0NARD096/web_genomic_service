from django.db import models
from django.contrib.auth.models import AbstractUser

class Profile(AbstractUser):
    email = models.EmailField(max_length=200, unique=True)
    phoneNumber = models.CharField(max_length=15)
    first_name = models.CharField(max_length=200)
    last_name = models.CharField(max_length=200)
    role = models.CharField(max_length=50)

    def __str__(self):
        return str(self.username)
    


class Genome(models.Model):
    sequence = models.TextField()
    species = models.CharField(max_length=255)
    description = models.CharField(max_length=255, default='Default Description')
    type = models.CharField(max_length=255, default='Default Type')

    def __str__(self):
        return self.species + " " + self.type 
     
""" 
class Annotation(models.Model):
    text = models.TextField()
    user = models.ForeignKey(User, on_delete=models.CASCADE, related_name='annotations')
    sequence = models.OneToOneField('Sequence', on_delete=models.SET_NULL, null=True, blank=True)
    isValidated = models.BooleanField(default=False)
    createdAt = models.DateTimeField(auto_now_add=True)
    
    def __str__(self):
        return f"Annotation by {self.user} - {self.text}"

class Sequence(models.Model):
    isUnique = models.BooleanField()
    genome = models.ForeignKey(Genome, on_delete=models.CASCADE, related_name='sequences')
    cdsFile = models.FileField(upload_to='cds_files/')
    pepFile = models.FileField(upload_to='pep_files/')
    
    def __str__(self):
        return f"Sequence for {self.genome}" """
    

from django.db import models
from django.contrib.auth.models import AbstractUser

class Profile(AbstractUser):
    email = models.EmailField('email',max_length=200, unique=True)
    phoneNumber = models.CharField('phone_number',max_length=15)
    first_name = models.CharField('first_name',max_length=200)
    last_name = models.CharField('last_name',max_length=200)
    role = models.CharField('role',max_length=50)
    is_approved = models.BooleanField('approved',default=False)

    def __str__(self):
        return str(self.username)
    
class AnnotationGenome(models.Model):
    species = models.CharField(max_length=255)
    upload_time = models.DateTimeField()
    annotated = models.BooleanField('annotated',default=False)
    annotator = models.ForeignKey(Profile, on_delete=models.CASCADE, null=True, blank=True)

    def __str__(self):
        return self.species

class Genome(models.Model):
    annotations = models.OneToOneField(
        AnnotationGenome, 
        on_delete=models.CASCADE
        )
    sequence = models.TextField()
    chromosome = models.CharField(max_length=255)
    start = models.IntegerField(null=True)
    end = models.IntegerField(null=True)
    
    def __str__(self):
        return self.chromosome 
    
class GeneProtein(models.Model):
    accession_number = models.CharField(max_length=255)
    sequence = models.TextField()
    type = models.CharField(max_length=255)
    start = models.IntegerField()
    end = models.IntegerField()
    genome = models.ForeignKey(Genome, on_delete=models.CASCADE)
    annotated = models.BooleanField('annotated',default=False)

    def __str__(self):
        return self.accession_number + " " + self.type 

""" 
class Annotation(models.Model):
    text = models.TextField()
    user = models.ForeignKey(User, on_delete=models.CASCADE, related_name='annotations')
    sequence = models.OneToOneField('Sequence', on_delete=models.SET_NULL, null=True, blank=True)
    isValidated = models.BooleanField(default=False)
    createdAt = models.DateTimeField(auto_now_add=True)
    
    def __str__(self):
        return f"Annotation by {self.user} - {self.text}"
"""

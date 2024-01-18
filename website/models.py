from django.db import models

class User(models.Model):
    email = models.EmailField(max_length = 200)
    password = models.CharField(max_length = 200)
    firstName = models.CharField(max_length = 200)
    lastName = models.CharField(max_length = 200)
    phoneNumber = models.IntegerField()
    
    def __str__(self):
        return self.firstName + " " + self.lastName

class Genome(models.Model):
    sequence = models.TextField()
    species = models.CharField(max_length=255)
    description = models.CharField(max_length=255)
    type = models.CharField(max_length=255)

    def __str__(self):
        return self.species + " " + self.type 
   
    def delete_everything(self):
        Genome.objects.all().delete()
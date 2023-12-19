from django.db import models

# Create your models here.

class SequenceFile(models.Model):
    file = models.FileField(upload_to='fa_files/')

    def __str__(self):
        return self.file.name
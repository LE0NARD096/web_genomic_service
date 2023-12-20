from django.db import models

# Create your models here.

class FastaFile(models.Model):
    file = models.FileField(upload_to='fasta_files/')

    def __str__(self):
        return self.file.name
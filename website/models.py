from django.db import models
from django.contrib.auth.models import AbstractUser
from django.contrib.contenttypes.models import ContentType
from django.contrib.contenttypes.fields import GenericForeignKey, GenericRelation
from phonenumber_field.modelfields import PhoneNumberField
from django.utils.text import slugify
from tinymce.models import HTMLField

class Profile(AbstractUser):
    USER = 'user'
    ANNOTATOR = 'annotator'
    VALIDATOR = 'validator'
    ADMIN = 'admin'

    ROLE_CHOICES = [
        (USER, 'user'),
        (ANNOTATOR, 'annotator'),
        (VALIDATOR, 'validator'),
        (ADMIN, 'admin'),
    ]

    USERNAME_FIELD = 'email'
    REQUIRED_FIELDS = ['username']

    email = models.EmailField('email', max_length=200, unique=True)
    phoneNumber = PhoneNumberField('phonenumber',null=True, blank=True, unique=True)
    first_name = models.CharField('first_name', max_length=200)
    last_name = models.CharField('last_name', max_length=200)
    role = models.CharField('role', choices=ROLE_CHOICES, max_length=50)
    is_approved = models.BooleanField('approved', default=False)

    def __str__(self):
        return str(self.username)

class AnnotationStatus(models.Model):
    VALIDATED = 'validated'
    PENDING = 'pending'
    REFUSED = 'refused'

    ROLE_CHOICES = [
        (VALIDATED, 'validated'),
        (PENDING, 'pending'),
        (REFUSED, 'refused'),
    ]
    content_type = models.ForeignKey(ContentType, 
                                     on_delete=models.CASCADE)
    object_id = models.PositiveIntegerField()
    validator = models.ForeignKey(Profile, 
                                  on_delete=models.CASCADE)
    content_object = GenericForeignKey('content_type', 'object_id')
    comment = models.TextField(blank=True)
    status = models.CharField('status', choices=ROLE_CHOICES, max_length=50)
    validation_time = models.DateTimeField(auto_now=True)

    def __str__(self):
        return str(self.object_id) + ': ' + self.status

    class Meta:
        verbose_name_plural = "Annotation status"


class Genome(models.Model):
    sequence = models.TextField(null=True,blank=True)
    chromosome = models.CharField(max_length=255)
    start = models.IntegerField(null=True)
    end = models.IntegerField(null=True)
    upload_time = models.DateTimeField(auto_now=True)
    is_validated = models.BooleanField('validated',default=False)
    
    def __str__(self):
        return self.chromosome 
    
class AnnotationGenome(models.Model):
    species = models.CharField(max_length=255, default='unknown')
    annotation_time = models.DateTimeField(null=True)
    annotator = models.ForeignKey(Profile, on_delete=models.CASCADE, null=True, blank=True)
    genome = models.OneToOneField(Genome, on_delete=models.CASCADE)
    is_annotated = models.BooleanField('annotated',default=False)
    comments = GenericRelation(AnnotationStatus)
    
    def __str__(self):
        return self.species

class GeneProtein(models.Model):
    accession_number = models.CharField(max_length=255)
    sequence = models.TextField()
    type = models.CharField(max_length=255)
    start = models.IntegerField()
    end = models.IntegerField()
    genome = models.ForeignKey(Genome, 
                               on_delete=models.CASCADE)
    upload_time = models.DateTimeField(auto_now=True)
    is_validated = models.BooleanField('validated',default=False)
   

    def __str__(self):
        return self.accession_number + " " + self.type 


class AnnotationProtein(models.Model):
    gene = models.CharField(max_length=255,default="unknown")
    transcript = models.CharField(max_length=255,default="unknown")
    gene_biotype = models.CharField(max_length=255,default="unknown")
    transcript_biotype = models.CharField(max_length=255,default="unknown")
    gene_symbol = models.CharField(max_length=255,default="unknown")
    description = models.CharField(max_length=1000,default="empty")
    annotation_time = models.DateTimeField(null=True)
    geneprotein = models.OneToOneField(GeneProtein, 
                                       on_delete=models.CASCADE, 
                                       null=True, 
                                       blank=True)
    annotator = models.ForeignKey(Profile, 
                                  on_delete=models.CASCADE, 
                                  null=True, 
                                  blank=True)
    is_annotated = models.BooleanField('annotated',default=False)
    comments = GenericRelation(AnnotationStatus)

    def __str__(self):
        return self.gene
comments = GenericRelation(AnnotationStatus)

class Post(models.Model):
    title = models.CharField(max_length = 500)
    url = models.SlugField(max_length=500, unique=True, blank=True)
    author = models.ForeignKey(Profile, 
                               on_delete=models.CASCADE)
    content = HTMLField()
    publication_date = models.DateTimeField(auto_now_add = True)
    comments = GenericRelation(AnnotationStatus)
    
    def __str__(self):
        return self.title
    
    def save(self, *args, **kwargs):
        if not self.url:
            self.url = slugify(self.title)
        super(Post, self).save(*args, **kwargs)



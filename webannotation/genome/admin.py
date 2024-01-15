from django.contrib import admin
from .models import User, Role, Annotation, Sequence, Genome, File

# Register your models here.
admin.site.register(Role)
admin.site.register(Annotation)
admin.site.register(Sequence)
admin.site.register(Genome)
admin.site.register(File)


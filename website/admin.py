from django.contrib import admin
from .models import Profile, Genome, GeneProtein

@admin.register(Profile)
class ProfileAdmin(admin.ModelAdmin):
    fields = (('username','is_approved'),'email','role')
    list_display = ('username','role','is_approved')
    ordering = ('-is_approved',)

@admin.register(Genome)
class GenomeAdmin(admin.ModelAdmin):
    fields = (('chromosome','annotated'),'start','end')
    list_display = ('chromosome','annotated')
    ordering = ('-annotated',)

@admin.register(GeneProtein)
class GenomeAdmin(admin.ModelAdmin):
    fields = (('accession_number','type'),'start','end')
    list_display = ('accession_number','genome','annotated')
    ordering = ('-annotated',)
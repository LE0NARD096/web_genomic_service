from django.contrib import admin
from .models import Profile, Genome, GeneProtein, AnnotationGenome
from django.utils.html import format_html

@admin.register(Profile)
class ProfileAdmin(admin.ModelAdmin):
    fields = (('username','is_approved'),'email','role')
    list_display = ('username','role','is_approved')
    ordering = ('-is_approved',)

@admin.register(GeneProtein)
class GenomeAdmin(admin.ModelAdmin):
    fields = (('accession_number','type'),'start','end','genome')
    list_display = ('accession_number','genome','annotated')
    ordering = ('-annotated',)

@admin.register(AnnotationGenome)
class AnnotationGenomeAdmin(admin.ModelAdmin):
    fields = (('species','annotator','annotated'),'upload_time')
    list_display = ('species','genome','annotated')
    ordering = ('-annotated',)

@admin.register(Genome)
class GenomeAdmin(admin.ModelAdmin):
    fields = (('chromosome', 'annotations'), 'start', 'end')
    list_display = ('chromosome', 'annotations_info', 'start', 'end')

    def annotations_info(self, obj):
        return obj.annotations.annotated
        
    annotations_info.boolean = True
    annotations_info.short_description = 'Annotated'
    
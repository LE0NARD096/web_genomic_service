from django.contrib import admin
from .models import Profile, Genome, GeneProtein, AnnotationGenome, AnnotationProtein

@admin.register(Profile)
class ProfileAdmin(admin.ModelAdmin):
    fields = (('username','is_approved'),'email','role')
    list_display = ('username','role','is_approved')
    ordering = ('-is_approved',)

@admin.register(GeneProtein)
class GeneProteinAdmin(admin.ModelAdmin):
    fields = (('accession_number','type','annotated'),'sequence','start','end','genome')
    list_display = ('accession_number','genome','annotated','upload_time')

@admin.register(Genome)
class GenomeAdmin(admin.ModelAdmin):
    fields = (('chromosome','annotated'), 'start', 'end')
    list_display = ('species', 'chromosome', 'annotated', 'upload_time')

    
    def species(self, obj):
        if obj.annotationgenome:
            return obj.annotationgenome.species
        return None

    species.short_description = 'species'  



@admin.register(AnnotationGenome)
class AnnotationGenomeAdmin(admin.ModelAdmin):
    fields = (('species','annotator','is_validated'),'annotation_time')
    list_display = ('species','genome','is_validated','annotator','annotation_time')

    actions = ['mark_as_annotated','mark_as_not_annotated']

    def mark_as_annotated(self, request, queryset):
        queryset.update(annotated=True)
        self.message_user(request, f'Marked {queryset.count()} items as annotated.')

    def mark_as_not_annotated(self, request, queryset):
        queryset.update(annotated=False)
        self.message_user(request, f"Marked {queryset.count()} items aren't annotated anymore.")


@admin.register(AnnotationProtein)
class AnnotationProteinAdmin(admin.ModelAdmin):
    list_display = ('gene', 'is_validated','annotator', 'annotation_time')
    actions = ['mark_as_approved']

    def mark_as_approved(self, request, queryset):
        queryset.update(annotated=True)
        self.message_user(request, f'Marked {queryset.count()} items as annotated.')



from django.contrib import admin
from .models import Profile, Genome, GeneProtein, AnnotationGenome, AnnotationProtein

@admin.register(Profile)
class ProfileAdmin(admin.ModelAdmin):
    fields = (('username','is_approved'),'email','role')
    list_display = ('username','role','is_approved')
    ordering = ('-is_approved',)

@admin.register(Genome)
class GenomeAdmin(admin.ModelAdmin):
    fields = ('chromosome', 'start', 'end')
    list_display = ('species', 'chromosome', 'annotated', 'upload_time')

    
    def species(self, obj):
        if obj.annotationgenome:
            return obj.annotationgenome.species
        return None

    species.short_description = 'species'  



@admin.register(AnnotationGenome)
class AnnotationGenomeAdmin(admin.ModelAdmin):
    fields = (('species','annotator','annotated'),'annotation_time')
    list_display = ('species','genome','annotations_info','annotator','annotation_time')

    def annotations_info(self, obj):
        return obj.genome.annotated
    
    annotations_info.boolean = True
    annotations_info.short_description = 'Annotated'
    

    actions = ['mark_as_annotated','mark_as_not_annotated']

    def mark_as_annotated(self, request, queryset):
        queryset.update(annotated=True)
        self.message_user(request, f'Marked {queryset.count()} items as annotated.')

    def mark_as_not_annotated(self, request, queryset):
        queryset.update(annotated=False)
        self.message_user(request, f"Marked {queryset.count()} items aren't annotated anymore.")



    

@admin.register(AnnotationProtein)
class AnnotationProteinAdmin(admin.ModelAdmin):
    list_display = ('gene', 'annotations_info', 'annotation_time')
    actions = ['mark_as_approved']

    def annotations_info(self, obj):
        return obj.geneprotein.annotated

    annotations_info.boolean = True
    annotations_info.short_description = 'Annotated'

    def mark_as_approved(self, request, queryset):
        queryset.update(annotated=True)
        self.message_user(request, f'Marked {queryset.count()} items as annotated.')



from .models import Genome, GeneProtein, AnnotationGenome, AnnotationProtein, Profile
from django.db.models import Q

def notifications(request):
    notification_validator = 0
    notification_annotator = 0

    if request.user.is_authenticated:
        if request.user.role == "validator" or request.user.is_superuser:
            notification_validator = AnnotationGenome.objects.select_related('genome').filter(is_annotated = True, genome__is_validated=False).count() \
                                    + AnnotationProtein.objects.select_related('geneprotein').filter(is_annotated = True, geneprotein__is_validated=False).count() \
                                    + Genome.objects.filter(Q(is_validated=False) &
                                                            Q(annotationgenome__isnull=True)).defer('sequence').count() \
                                    + GeneProtein.objects.filter(Q(is_validated=False) &
                                                                 Q(annotationprotein__isnull=True)).defer('sequence').count()
                                    
        if request.user.role == "annotator" or request.user.is_superuser:
            annotator = Profile.objects.get(pk=request.user.id)
            notification_annotator = AnnotationGenome.objects.filter(is_annotated=False, annotator=annotator,
                                                                    genome__sequence__isnull=False).count() \
                                     + AnnotationProtein.objects.filter(is_annotated=False, annotator=annotator).count()

    return {
        'notifications_validator': notification_validator,
        'notifications_annotator': notification_annotator
    }

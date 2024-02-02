from .models import Genome, GeneProtein, AnnotationGenome, AnnotationProtein, Profile

def notifications(request):
    notification_validator = 0
    notification_annotator = 0

    if request.user.is_authenticated:
        if request.user.role == "validator" or request.user.is_superuser:
            notification_validator = Genome.objects.filter(is_validated=False, annotationgenome__isnull=True).count() \
                                    + GeneProtein.objects.filter(is_validated=False, annotationprotein__isnull=True).count()

        if request.user.role == "annotator" or request.user.is_superuser:
            annotator = Profile.objects.get(pk=request.user.id)
            notification_annotator = AnnotationGenome.objects.filter(is_annotated=False, annotator=annotator,
                                                                    genome__sequence__isnull=False).count() \
                                     + AnnotationProtein.objects.filter(is_annotated=False, annotator=annotator).count()

    return {
        'notifications_validator': notification_validator,
        'notifications_annotator': notification_annotator
    }

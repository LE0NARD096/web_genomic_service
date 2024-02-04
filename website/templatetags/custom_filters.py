from django import template

register = template.Library()

@register.filter(name='can_edit_annotation')
def can_edit_annotation(annotator, request_user):
    return annotator.id == request_user.id and (request_user.role == "annotator" or request_user.is_superuser)

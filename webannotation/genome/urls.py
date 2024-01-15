from django.urls import path, include
from . import views

urlpatterns = [
    path('', views.home, name='home'),
    path('login/', views.login_view, name='login_url'),
    path('create_annotation/', views.create_annotation, name='create_annotation'),
    path('edit_annotation/<int:annotation_id>/', views.edit_annotation, name='edit_annotation'),
    path('validate_annotation/<int:annotation_id>/', views.validate_annotation, name='validate_annotation'),
    path('list_sequences/', views.list_sequences, name='list_sequences'),
    path('view_annotation/<int:annotation_id>/', views.view_annotation, name='view_annotation'),
    
]

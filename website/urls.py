from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name='home'),
    path('login/', views.login_view, name='login'), 
    path('logout/', views.logout_view, name='logout'), 
    path('search/', views.search_results, name='search_results'),
    path('upload/', views.upload, name='upload_file'),
    path('register/', views.register_view, name='register'),
    path('text_extraction/', views.text_extraction, name='text_extraction'),
    path('visualisation_<str:type>/<int:id>/', views.visualisation_sequence, name='visualisation'),
    path('validator_view/', views.validator_view,name='validator_view'),
    path('validate_annotators/',views.assigned_annotators,name='validate_annotators'),
    path('annotator_dashboard/', views.annotator_view, name='annotator_dashboard'),
    path('annotation_<str:type_of_sequence>/<int:id>/', views.sequence_view, name='annotation'),
    path('validation_<str:type_of_sequence>/<int:id>/', views.sequence_view, name='validation'),
    path('annotated/',views.save_annotation,name='annotated'),
    path('validated/',views.validate_include_database,name='validated'),
    path('blast_results/', views.blast_search, name='blast_results'),
    path('delete_<str:type_of_sequence>/<int:id>/', views.delete_sequence_from_database, name='delete'),
    path('include_sequence_database_<str:sequence_type>/<int:id_get_sequence>/',views.validate_include_database,name='include'),
    path('visual_annotation_<str:sequence_type>/<int:selected_genome_id>/', views.visualize_genome_genes, name = "visual_annotation"),
]



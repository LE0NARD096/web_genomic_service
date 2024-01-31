from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name='home'),
    path('search/', views.search_results, name='search_results'),
    path('upload/', views.upload, name='upload_file'),
    path('login/', views.login_view, name='login'),
    path('logout.html/', views.logout_view, name='logout'),  
    path('register/', views.register_view, name='register'),
    path('text_extraction/', views.text_extraction, name='text_extraction'),
    path('visualisation_<str:type>/<int:id>/', views.visualisation_sequence, name='visualisation'),
    path('validator_view', views.validator_view,name='validator_view'),
    path('validate_annotators',views.assigned_annotators,name='validate_annotators'),
    path('annotator_dashboard', views.annotator_view, name='annotator_dashboard'),
    path('annotation_<str:type>/<int:id>/', views.annotator_modify, name='annotation'),
    path('annotated',views.save_annotation,name='annotated')
]

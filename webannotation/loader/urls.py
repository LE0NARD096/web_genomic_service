from django.urls import path
from .views import fasta_file_view

urlpatterns = [
    path('fasta/', fasta_file_view, name='fasta_file_view'),
]

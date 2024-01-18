from django import forms
from django.forms import ModelForm
from .models import Genome

class GenomeSearchForm(forms.Form):
    sequence = forms.CharField(label='Sequence query', required=True,widget=forms.Textarea)
    species = forms.CharField(label='Species', required=False)
    output_type = forms.ChoiceField(label='Search in', choices=[('genome', 'Génome'), ('gene_protein', 'Gène/Protéine')])

class Upload_data(forms.Form):
    sequence = forms.FileField(help_text="Upload a fasta file",allow_empty_file=False)
    species = forms.CharField(max_length=200)


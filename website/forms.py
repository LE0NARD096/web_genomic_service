from typing import Any
from django import forms
from .models import Profile, AnnotationProtein, AnnotationGenome
from django.contrib.auth.forms import UserCreationForm
from django.core.exceptions import ValidationError
from Bio import SeqIO
from io import StringIO

class GenomeSearchForm(forms.Form):
    sequence = forms.CharField(label='Sequence query', required=True,widget=forms.Textarea)
    species =forms.CharField(label='Species', required=False)
    chromosome = forms.CharField(label='Chromosome', required=False)
    transcript = forms.CharField(label='Transcript', required=False)
    gene = forms.CharField(label='Gene', required=False)
    output_type = forms.ChoiceField(label='Search in', choices=[('genome', 'Génome'), ('gene_protein', 'Gène/Protéine')])
    
    def clean(self):
        """
        Validate the form input
        """
        cleaned_data = super().clean()
        type = cleaned_data.get('output_type')

        #if type == 'gene_protein':
            #gene = cleaned_data.get('gene')
           # if not gene:
                #self.add_error('gene', ("This field is required for gene/protein"))
    

        return cleaned_data


class GenomeAnnotate(forms.ModelForm):
    class Meta:
        model = AnnotationGenome
        fields = '__all__'

class DownloadTextForm(forms.Form):
    StartPosition = forms.IntegerField(label='Start', required=False)
    EndPosition = forms.IntegerField(label='End', required=False)
    species = forms.CharField(label='Species', required=False)
    output_type = forms.ChoiceField(label='Search in', choices=[('genome', 'Génome'), ('gene_protein', 'Gène/Protéine')])


class Upload_data(forms.Form):
    sequence = forms.FileField(
        help_text="Upload a fasta file",
        allow_empty_file=False,
        widget=forms.ClearableFileInput(attrs={'accept': '.fa'})
    )
    output_type = forms.ChoiceField(label='Upload in', choices=[('genome', 'Genome'), ('gene_protein', 'Gene/Protein')])
    annotated = forms.BooleanField(label='Annotated sequence', required=False)



class UserRegistrationForm(UserCreationForm):
    email = forms.EmailField(label='Email', required=True)
    first_name = forms.CharField(label='First Name', required=True)
    last_name = forms.CharField(label='Last Name', required=True)
    phone_number = forms.CharField(label='Phone Number', required=False)

    ROLE_CHOICES = [
        ('user', 'User'),
        ('annotator', 'Annotator'),
        ('validator', 'Validator'),
    ]

    role = forms.ChoiceField(label='Role', choices=ROLE_CHOICES)

    class Meta:
        model = Profile
        fields = ['username','email', 'first_name', 'last_name', 'phone_number', 'role', 'password1', 'password2']

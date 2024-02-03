from typing import Any
from django import forms
from .models import Profile, AnnotationProtein, AnnotationGenome, GeneProtein, Genome
from django.contrib.auth.forms import UserCreationForm
from django.core.exceptions import ValidationError
from Bio import SeqIO
from io import StringIO
from django.utils.safestring import mark_safe
from django.urls import reverse



class GenomeSearchForm(forms.Form):
    sequence = forms.CharField(label='Sequence:', required=True,widget=forms.Textarea)
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

class ProteinAnnotate(forms.ModelForm):
    class Meta:
        model = AnnotationProtein
        fields = '__all__'
        exclude = ['annotator','is_annotated','annotation_time']

class SequenceProtein(forms.ModelForm):
    class Meta:
        model = GeneProtein
        fields = '__all__'
        exclude = ['is_validated']

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        for field_name in ['sequence']:  
            self.fields[field_name].label = ''

class GenomeAnnotate(forms.ModelForm):
    class Meta:
        model = AnnotationGenome
        fields = '__all__'
        exclude = ['annotator','is_annotated','annotation_time','genome']

class SequenceGenome(forms.ModelForm):
    class Meta:
        model = Genome
        fields = '__all__'
        exclude = ['is_validated','sequence']



class DownloadTextForm(forms.Form):
    output_type = forms.ChoiceField(label='Search in', choices=[('genome', 'Génome'), ('gene_protein', 'Gène/Protéine')])
    start_position = forms.IntegerField(label='Start', required=False)
    end_position = forms.IntegerField(label='End', required=False)
    chromosome = forms.CharField(label='Chromosome', required=False)
    species = forms.CharField(label='Species', required=False)
    ##
    gene = forms.CharField(label='Gene', required=False)
    transcript = forms.CharField(label='Transcript', required=False)
    function = forms.CharField(label='Function', required=False)


class Upload_data(forms.Form):
    sequence = forms.FileField(
        help_text="Upload a fasta file",
        allow_empty_file=False,
        widget=forms.ClearableFileInput(attrs={'class': 'form-control', 'placeholder': 'Upload a fasta file', 'accept': '.fa'})
    )
    output_type = forms.ChoiceField(label='Upload in', 
                                    choices=[('genome', 'Genome'), ('gene_protein', 'Gene/Protein')],
                                    widget=forms.Select(attrs={'class':'regDropDown'}))
    annotated = forms.BooleanField(required=False,
                                   widget=forms.CheckboxInput(attrs={'class': 'form-check-input'}))



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

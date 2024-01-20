from django import forms
from .models import Profile
from django.contrib.auth.forms import UserCreationForm




class GenomeSearchForm(forms.Form):
    sequence = forms.CharField(label='Sequence query', required=True,widget=forms.Textarea)
    species = forms.CharField(label='Species', required=False)
    output_type = forms.ChoiceField(label='Search in', choices=[('genome', 'Génome'), ('gene_protein', 'Gène/Protéine')])

class Upload_data(forms.Form):
    sequence = forms.FileField(help_text="Upload a fasta file",allow_empty_file=False)
    species = forms.CharField(max_length=200)


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

""" class GenomeForm(forms.ModelForm):
    class Meta:
        model = Genome
        fields = ['sequence', 'species', 'description', 'type'] """

""" class AnnotationForm(forms.ModelForm):
    class Meta:
        model = Annotation
        fields = ['text', 'sequence', 'isValidated']


class SequenceForm(forms.ModelForm):
    class Meta:
        model = Sequence
        fields = ['isUnique', 'genome', 'cdsFile', 'pepFile'] """
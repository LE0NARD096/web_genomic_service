from django import forms
from django.forms import ModelForm
from .models import Genome, Annotation, Sequence, User
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
    phone_number = forms.CharField(label='Phone Number', required=True)
    first_name = forms.CharField(label='First Name', required=True)
    last_name = forms.CharField(label='Last Name', required=True)
    username = forms.CharField(label='Username', required=True)  # Add this line

    ROLE_CHOICES = [
        ('user', 'User'),
        ('annotator', 'Annotator'),
        ('validator', 'Validator'),
    ]

    role = forms.ChoiceField(label='Role', choices=ROLE_CHOICES)

    class Meta(UserCreationForm.Meta):
        model = User
        fields = ['username', 'email', 'password1', 'password2', 'first_name', 'last_name', 'phone_number', 'role']
        widgets = {
            'password1': forms.PasswordInput(),
            'password2': forms.PasswordInput(),
        }

class LoginForm(forms.Form):
    username = forms.CharField(max_length=65)
    password = forms.CharField(max_length=65, widget=forms.PasswordInput)

class GenomeForm(forms.ModelForm):
    class Meta:
        model = Genome
        fields = ['sequence', 'species', 'description', 'type']

class AnnotationForm(forms.ModelForm):
    class Meta:
        model = Annotation
        fields = ['text', 'sequence', 'isValidated']


class SequenceForm(forms.ModelForm):
    class Meta:
        model = Sequence
        fields = ['isUnique', 'genome', 'cdsFile', 'pepFile']

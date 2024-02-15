from typing import Any
from django import forms
from .models import Profile, AnnotationProtein, AnnotationGenome, GeneProtein, Genome, AnnotationStatus, Post
from django.contrib.auth.forms import UserCreationForm
from django.core.exceptions import ValidationError
from Bio import SeqIO
from io import StringIO
from django.utils.safestring import mark_safe
from django.urls import reverse
import re
from tinymce.widgets import TinyMCE

class CommentForm(forms.Form):
    
    comment = forms.CharField(required=False,widget=forms.Textarea)
    
    STATUS_CHOICES = [
        ('validated', 'Validate'),
        ('refused', 'Refuse'),
    ]
    
    status = forms.ChoiceField(label='Status', choices=STATUS_CHOICES)

    def __init__(self, *args, **kwargs):
        super(CommentForm, self).__init__(*args, **kwargs)

        self.fields["comment"].label = ''
   

class UpdateForm(forms.ModelForm):
    class Meta:
        model = Profile
        fields = ['username','first_name', 'last_name', 'phoneNumber']
    

class GenomeSearchForm(forms.Form):
    output_type = forms.ChoiceField(label='Search in', choices=[('genome', 'Genome'), ('gene_protein', 'Gene/Protein')])
    sequence = forms.CharField(label='Sequence:', required=True,widget=forms.Textarea)
    species =forms.CharField(label='Species', required=False)
    chromosome = forms.CharField(label='Chromosome', required=False)
    ##
    transcript = forms.CharField(label='Transcript', required=False)
    gene = forms.CharField(label='Gene', required=False)
    function = forms.CharField(label='Function',required=False)
    database = forms.ChoiceField(label = 'Database', choices=[('BactaHub', 'BactaHub'), ('NCBI','NCBI')])
    
    def clean_sequence(self):
        sequence = self.cleaned_data.get('sequence')
        type_sequence = self.cleaned_data.get('output_type')

        if type_sequence == 'genome' and re.fullmatch(r'^[ATGCN]*$', sequence) is None:
           raise ValidationError(("You can't search a protein fragment in a genome"), code='invalid')

        if len(sequence) < 3:
           raise ValidationError(('Sequence too short'), code='invalid')
        
        if len(sequence) > 924:
            raise ValidationError(('Sequence too long'), code='invalid')

        if re.search('%',sequence):
            if re.fullmatch(r'%[ATGCN]+%[ATGCN]+%',sequence) is None:
                if re.fullmatch(r'%[ARNDCQEGHILKMFPSTWYV]+%[ARNDCQEGHILKMFPSTWYV]+%',sequence) is None:
                    raise ValidationError(("The protein sequence has a" + ' % ' + "symbol but doesn't respect the pattern " + "%" + "any protein/DNA" + "%" + "any protein/DNA" + "%"), code='invalid')
        
        else:
            if re.fullmatch(r'^[ATGCN]*$',sequence) is None:
                if re.fullmatch(r'^[ACDEFGHIKLMNPQRSTVWY]*$',sequence) is None:
                    raise ValidationError(("The sequence does not respect a protein or DNA format"), code='invalid')
            
        return sequence
        

class ProteinAnnotate(forms.ModelForm):
    class Meta:
        model = AnnotationProtein
        fields = '__all__'
        exclude = ['annotator','is_annotated','annotation_time', 'geneprotein']
        read_only_fields = ('is_active', 'is_staff')
    
    def __init__(self, *args, **kwargs):
        self.current_user = kwargs.pop('current_user', None)
        super(ProteinAnnotate, self).__init__(*args, **kwargs)

        if self.current_user.role == 'validator':
            
            for field_name, field in self.fields.items():
                field.widget.attrs['readonly'] = True

class SequenceProtein(forms.ModelForm):
    class Meta:
        model = GeneProtein
        fields = '__all__'
        exclude = ['is_validated']

    def __init__(self, *args, **kwargs):
        self.current_user = kwargs.pop('current_user', None)
        super(SequenceProtein, self).__init__(*args, **kwargs)

        self.fields["sequence"].label = ''

        if self.current_user.role == 'validator':
            self.fields['genome'].widget.attrs['disabled'] = True
            self.fields['genome'].required = False
            for field_name, field in self.fields.items():
                field.widget.attrs['readonly'] = True
        

class GenomeAnnotate(forms.ModelForm):
    class Meta:
        model = AnnotationGenome
        fields = '__all__'
        exclude = ['annotator','is_annotated','annotation_time','genome']
    
    def __init__(self, *args, **kwargs):
        self.current_user = kwargs.pop('current_user', None)
        super(GenomeAnnotate, self).__init__(*args, **kwargs)


        if self.current_user.role == 'validator':
            for field_name, field in self.fields.items():
                field.widget.attrs['readonly'] = True

class SequenceGenome(forms.ModelForm):
    class Meta:
        model = Genome
        fields = '__all__'
        exclude = ['is_validated','sequence']

    def __init__(self, *args, **kwargs):
        self.current_user = kwargs.pop('current_user', None)
        super(SequenceGenome, self).__init__(*args, **kwargs)

        if self.current_user.role == 'validator':
            for field_name, field in self.fields.items():
                field.widget.attrs['readonly'] = True
    

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



class ReplyForm(forms.Form):
    comment = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Post a comment'}))

class CreatePostForm(forms.ModelForm):
    class Meta:
        model = Post
        fields = ['title','content']
        widgets = {'content': TinyMCE(mce_attrs={'height': 280}), 'title': forms.TextInput(attrs={'placeholder': 'Title'})}

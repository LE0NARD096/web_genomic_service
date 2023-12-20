from django import forms
from .models import FastaFile

class FastaFileForm(forms.ModelForm):
    class Meta:
        model = FastaFile
        fields = ['file']

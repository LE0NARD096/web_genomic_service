from django import forms

class GenomeSearchForm(forms.Form):
    sequence = forms.CharField(label='Séquence', required=False)
    species = forms.CharField(label='Espèce', required=False)
    output_type = forms.ChoiceField(label='Type de sortie', choices=[('genome', 'Génome'), ('gene_protein', 'Gène/Protéine')])

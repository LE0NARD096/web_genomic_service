from django.shortcuts import render
from .models import User, Genome
from .forms import GenomeSearchForm
from Bio import SeqIO
from Bio.Seq import Seq


def home(request):
    all_users = User.objects.all
    return render(request,'home.html',{'all': all_users})

def search_results(request):
    if request.method == 'POST':
        form = GenomeSearchForm(request.POST)
        if form.is_valid():
            sequence_query = form.cleaned_data['sequence']
            species_query = form.cleaned_data['species']
            output_type = form.cleaned_data['output_type']

            if output_type == 'genome':
                    results = Genome.objects.filter(sequence__contains=sequence_query, species__contains=species_query)
           
            return render(request, 'search_results.html', {'results': results, 'output_type': output_type})
    else:
        form = GenomeSearchForm()

    return render(request, 'search_form.html', {'form': form})


from django.shortcuts import render,redirect
from .models import User, Genome
from .forms import GenomeSearchForm,Upload_data
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO

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

            print(request.POST)
            if output_type == 'genome':
                results = Genome.objects.filter(species__contains=species_query,type = 'genome')
            else:
                results = Genome.objects.filter(species__contains=species_query, type__in=['pep', 'cds'])
            c = 0
            final_result = []
            for DNAsequence in results:
                c += 1  # Use the correct += operator
                my_dna = Seq(DNAsequence.sequence)
                if my_dna.count(sequence_query) > 0:
                    final_result.append(DNAsequence)
            
            print(final_result)

            return render(request, 'search_results.html', {'results': final_result})
    else:
        form = GenomeSearchForm()

    return render(request, 'search_form.html', {'form': form})

def upload(request):
    if request.POST:
        form = Upload_data(request.POST, request.FILES)
        if form.is_valid():
            file = request.FILES['sequence']
            file_content = file.read().decode('utf-8')
            file_stream = StringIO(file_content)
            for record in SeqIO.parse(file_stream, 'fasta'):
                    sequence = record.seq
                    if record.description.split()[0].lower() == 'chromosome':
                         type = 'genome'
                    else:
                         type = record.description.split()[1]
                    
                    species = request.POST['species']
                    description = record.description
                    genome = Genome.objects.create(sequence=sequence,species=species,description=description,type=type)
                    genome.save()
            return redirect(home)
    return render(request,'upload.html',{'form':Upload_data})


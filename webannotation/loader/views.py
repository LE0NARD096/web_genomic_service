from django.shortcuts import render, redirect
from .models import FastaFile
from .forms import FastaFileForm
from Bio import SeqIO

def fasta_file_view(request):
    if request.method == 'POST':
        form = FastaFileForm(request.POST, request.FILES)
        if form.is_valid():
            form.save()
            return redirect('fasta_file_view')

    else:
        form = FastaFileForm()

    fasta_files = FastaFile.objects.all()
    sequences = []

    for fasta_file in fasta_files:
        file_path = fasta_file.file.path
        with open(file_path) as fasta:
            sequences.extend(SeqIO.parse(fasta, "fasta"))

    return render(request, 'fasta_viewer.html', {'form': form, 'sequences': sequences})

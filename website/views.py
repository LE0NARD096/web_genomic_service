from django.shortcuts import render,redirect, get_object_or_404
from .models import Profile, Genome, GeneProtein, AnnotationGenome, AnnotationProtein
from .forms import GenomeSearchForm,Upload_data, DownloadTextForm, GenomeAnnotate
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO
from django.contrib.auth import authenticate, login, logout
from django.contrib import messages
from django.http import HttpResponse
from Bio.SeqFeature import SeqFeature, FeatureLocation, SimpleLocation
from django.urls import reverse
import re
from django.contrib.auth.forms import AuthenticationForm, UserCreationForm
from django.contrib.auth.decorators import login_required
from django.shortcuts import redirect
from .forms import UserRegistrationForm
from django.utils import timezone
from django.db.models import Q

def home(request):
    register_url = reverse('register') 
    login_url = reverse('login')  

    return render(request, 'home.html', {'register': register_url, 'login': login_url})

def text_extraction(request):
    if request.method == 'POST':
        form = DownloadTextForm(request.POST)
        if form.is_valid():
            startposition = form.cleaned_data['StartPosition']
            endposition = form.cleaned_data['EndPosition']
            species_query = form.cleaned_data['species']
            output_type = form.cleaned_data['output_type']

            if output_type == 'genome':
                results = Genome.objects.filter(species__contains=species_query,type = 'genome')
            else:
                results = Genome.objects.filter(species__contains=species_query, type__in=['pep', 'cds'])
        
            lines = []
            flag_None = False

            if startposition == None:
                startposition = 0
            
            if endposition == None:
                flag_None = True
            
            if not results.exists():
                lines.append(f'No results')
                
            else:
                for type_sequence in results:
                    my_sequence = Seq(type_sequence.sequence)

                    if flag_None == True:
                        endposition = len(my_sequence)
                

                    print(startposition,endposition)
                    if type_sequence.type == 'pep':
                        f = SeqFeature(FeatureLocation(startposition, endposition), type="domain")
                        if f:
                            extraction = f.extract(my_sequence)
                            lines.append(f'>{type_sequence.description}\n{extraction}\n')
                    else :
                        f = SeqFeature(FeatureLocation(startposition, endposition), type="CDS")
                        if f:
                            extraction = f.extract(my_sequence)
                            lines.append(f'>{type_sequence.description}\n{extraction}\n')

            response = HttpResponse(content_type='text/plain')
            response['Content-Disposition'] = f'attachment; filename= {species_query}_{output_type}_query.txt'
            response.writelines(lines)
            return response
    else:
        form = DownloadTextForm()

    return render(request, 'text_extraction.html', {'form': form})


def register_view(request):
    if request.method == 'POST':
        form = UserRegistrationForm(request.POST)
        print(form)
        if form.is_valid():
            user = form.save(commit=False)
            user.username = user.username.lower()
            user.save()
            return render(request,'home.html')
    else:
        messages.error(request, 'Error creating your account or invalid login credentials. Please correct the errors below.')
        form = UserRegistrationForm()

    return render(request, 'Authentication/register.html', {'form': form})


def login_view(request):
    if request.method == 'POST':
        form = AuthenticationForm(request, request.POST)
        if form.is_valid():
            username = form.cleaned_data['username']
            password = form.cleaned_data['password']
            user = authenticate(request, username=username, password=password,is_approved=True)
           
            if user is not None:
                if user.is_approved:
                    login(request, user)
                    return redirect('/')
                else:
                     messages.error(request, 'Your account is not approved yet. Please wait for approval.')

    else:
        form = AuthenticationForm()
    
    return render(request, 'Authentication/login.html', {'form': form})

def logout_view(request):
    if request.method == 'POST':
        logout(request)
        messages.success(request, f'You have been logged out.')
        return redirect('/')
    return render(request,'Authentication/logout.html',{})

def search_results(request):
    if request.method == 'POST':
        form = GenomeSearchForm(request.POST)
        if form.is_valid():
            sequence_query = form.cleaned_data['sequence']
            species = form.cleaned_data['species']
            gene = form.cleaned_data['gene']
            transcript = form.cleaned_data['transcript']
            chromosome = form.cleaned_data['chromosome']
            output_type = form.cleaned_data['output_type']

           
           
            c = 0
            final_result = []

            if output_type == 'genome':
                results = Genome.objects.filter(
                            Q(annotated=True) &
                            Q(chromosome=chromosome) &
                            Q(annotationgenome__species__contains=species)
                            )

            else:
                 if validate(sequence_query,'dna'):
                    results = GeneProtein.objects.filter(
                            Q(annotated=True) &
                            Q(type='cds') &
                            Q(annotationprotein__gene__contains=gene) &
                            Q(annotationprotein__transcript__contains=transcript) &
                            Q(genome__chromosome__contains=chromosome) &
                            Q(genome__annotationgenome__species__contains=species)
                            )

                 else:
                     results = GeneProtein.objects.filter(
                            Q(annotated=True) &
                            Q(type='pep') &
                            Q(annotationprotein__gene__contains=gene) &
                            Q(annotationprotein__transcript__contains=transcript) &
                            Q(genome__chromosome__contains=chromosome) &
                            Q(genome__annotationgenome__species__contains=species)
                            )
            

            print(results)
            for DNAsequence in results:
                c += 1 
                my_dna = Seq(DNAsequence.sequence)
                if my_dna.count(sequence_query) > 0:
                    if output_type == 'genome':
                        result_dic = {
                            'type': 'genome',
                            'chromosome': DNAsequence.chromosome,
                            'id': DNAsequence.id,
                            }
                    else:
                        result_dic = {
                            'type': DNAsequence.type,
                            'chromosome': DNAsequence.accession_number,
                            'id': DNAsequence.id,
                            }
                    final_result.append(result_dic)
            
            return render(request, 'search_results.html', {'results': final_result})
    else:
        form = GenomeSearchForm()

    return render(request, 'search_form.html', {'form': form})

def validate(seq, alphabet='dna'):
    
    alphabets = {'dna': re.compile('^[acgtn]*$', re.I), 
             'protein': re.compile('^[acdefghiklmnpqrstvwy]*$', re.I)}


    if alphabets[alphabet].search(seq) is not None:
         return True
    else:
         return False

def upload(request):
    if request.POST:
        form = Upload_data(request.POST, request.FILES)
        if form.is_valid():
            file = request.FILES['sequence']
            type = request.POST['output_type']
            file_content = file.read().decode('utf-8')
            file_stream = StringIO(file_content)

            if type == 'genome':
                record = SeqIO.parse(file_stream, 'fasta')
                genome = next(record, None)
                
                sequence = genome.seq
                description = genome.description.split(':')
                chromosome = description[2]
                start = description[4]
                end = description[5]

                if Genome.objects.filter(chromosome=chromosome).exists():
                    genome = Genome.objects.get(chromosome=chromosome)
                    genome.sequence = sequence
                    genome.start = start
                    genome.end = end
                    
                    
                else:
                    
                    genome = Genome.objects.create(sequence=sequence,
                                               chromosome=chromosome,
                                               start=start,
                                               end=end,
                                               upload_time=timezone.now())
                
                
                genome.save()

            else :
                for record in SeqIO.parse(file_stream, 'fasta'):
                    a_n = record.id
                    sequence = record.seq
                    type = record.description.split()[1].lower()
                    description = record.description.split(':')
                    start = description[3]
                    end = description[4]
                    chromosome = description[1]


                    if Genome.objects.filter(chromosome=chromosome).exists():
                        genome = Genome.objects.get(chromosome=chromosome)
                    else:
                        genome = Genome.objects.create(sequence='empty',
                                                       chromosome=chromosome,
                                                       upload_time=timezone.now())
                        genome.save()
                    
                    gene_protein = GeneProtein.objects.create(accession_number=a_n,
                                                        sequence=sequence,
                                                        type=type,
                                                        start=start,
                                                        end=end,
                                                        genome=genome,
                                                        upload_time = timezone.now())
                    gene_protein.save()

            return redirect(home)
        else:
            return render(request,'upload.html',{'form':form})

    return render(request,'upload.html',{'form':Upload_data})


def visualisation_sequence(request,type,id):
    if type == "genome":
        genome = Genome.objects.get(id=id)
    else:
        genome = GeneProtein.objects.get(id=id,type=type)
    return render(request, 'visualisation_sequence.html', {'genome_id': genome})


@login_required
def validator_view(request):
    unvalidated_genomes = Genome.objects.filter(annotated=False, annotationgenome__isnull=True)
    unvalidated_proteins = GeneProtein.objects.filter(annotated=False,annotationprotein__isnull=True)
    annotators = None


    if unvalidated_genomes != None or unvalidated_proteins != None:
        annotators = Profile.objects.filter(role='annotator')
    
        

    context = {
        'annotators': annotators,
        'unvalidated_genomes': unvalidated_genomes,
        'unvalidated_proteins': unvalidated_proteins
    }

    return render(request, 'Validator/validator_view.html', context)

@login_required
def assigned_annotators(request):
    if request.method == 'POST':
        print(request.POST)
        anno_id_protein = request.POST.getlist('annotation_id_protein')
        anno_id_genome = request.POST.getlist('annotation_id_genome')

        annotator_protein = request.POST.getlist('annotator_protein')
        annotator_genome = request.POST.getlist('annotator_genome')
        
        print(anno_id_genome)

        for i in range(len(anno_id_protein)):
            new_annotation_protein = AnnotationProtein.objects.create(
                annotator_id=annotator_protein[i],
                geneprotein_id=anno_id_protein[i]
            )
            new_annotation_protein.save()

        for i in range(len(anno_id_genome)):
            new_annotation_genome = AnnotationGenome.objects.create(
                annotator_id=annotator_genome[i],
                genome_id=anno_id_genome[i]
            )
            new_annotation_genome.save()
    
        return render(request, 'Validator/validate_annotators.html')

    return redirect('home')



@login_required
def annotator_view(request):
    annotator = request.user.id
    unvalidated_genome_annotations= Genome.objects.select_related('annotations').all()

    print(unvalidated_genome_annotations[1].annotations.annotator)
    
    AnnotationGenome.objects.filter(annotated=False, annotator=annotator)
    unvalidated_protein_annotations= AnnotationProtein.objects.filter(annotated=False, annotator=annotator)

    context = {
        'protein_annotations': unvalidated_protein_annotations,
        'genome_annotations': unvalidated_genome_annotations,
    }

    return render(request, 'Annotator/Annotator_dashboard.html', context)

@login_required
def annotator_modify(request,type,id):
    queryset = get_object_or_404(AnnotationGenome, id=id)
    form = GenomeAnnotate()
    form.fields = queryset
    return render(request,'Annotator/Annotate_sequences.html',{'form': GenomeAnnotate})








from django.shortcuts import render
from django.template.loader import render_to_string
from django.http import JsonResponse
from .models import GeneProtein
from Bio import SeqIO
import plotly.graph_objects as go

def visualizza_geni_genoma(genoma):
    # Recupera tutti i geni annotati sul genoma specificato
    geni_annotati = GeneProtein.objects.filter(annotated=True, genome=genoma)

    # Carica il genoma di riferimento utilizzando BioPython
    with open(genoma.genome_file.path, "r") as fasta_file:
        genoma_record = SeqIO.read(fasta_file, "fasta")

    # Crea una lista di tracce per i geni
    tracce = []
    for gene in geni_annotati:
        traccia_gene = go.Scatter(x=[gene.start, gene.end],
                                   y=[0, 0],
                                   mode='markers+text',
                                   marker=dict(size=10, symbol='line-ns-open', color='blue'),
                                   text=[gene.accession_number],
                                   name=gene.accession_number)
        tracce.append(traccia_gene)

    # Crea una traccia per il genoma di riferimento
    traccia_genoma = go.Scatter(x=[0, len(genoma_record.seq)],
                                y=[0, 0],
                                mode='lines',
                                line=dict(color='black', width=5),
                                name='Genoma di riferimento')

    # Aggiungi tutte le tracce al layout
    layout = go.Layout(title='Geni sul genoma',
                       xaxis=dict(title='Posizione nel genoma'),
                       yaxis=dict(visible=False),
                       showlegend=True)

    # Crea la figura
    figura = go.Figure(data=[traccia_genoma] + tracce, layout=layout)

    # Converti la figura in HTML
    grafico_html = figura.to_html(full_html=False)

    return grafico_html

def search_results(request):
    if request.method == 'POST':
        form = GenomeSearchForm(request.POST)
        if form.is_valid():
            sequence_query = form.cleaned_data['sequence']
            species = form.cleaned_data['species']
            gene = form.cleaned_data['gene']
            transcript = form.cleaned_data['transcript']
            chromosome = form.cleaned_data['chromosome']
            output_type = form.cleaned_data['output_type']

            c = 0
            final_result = []

            if output_type == 'genome':
                results = Genome.objects.filter(
                            Q(annotated=True) &
                            Q(chromosome=chromosome) &
                            Q(annotationgenome__species__contains=species)
                            )

            else:
                 if validate(sequence_query,'dna'):
                    results = GeneProtein.objects.filter(
                            Q(annotated=True) &
                            Q(type='cds') &
                            Q(annotationprotein__gene__contains=gene) &
                            Q(annotationprotein__transcript__contains=transcript) &
                            Q(genome__chromosome__contains=chromosome) &
                            Q(genome__annotationgenome__species__contains=species)
                            )

                 else:
                     results = GeneProtein.objects.filter(
                            Q(annotated=True) &
                            Q(type='pep') &
                            Q(annotationprotein__gene__contains=gene) &
                            Q(annotationprotein__transcript__contains=transcript) &
                            Q(genome__chromosome__contains=chromosome) &
                            Q(genome__annotationgenome__species__contains=species)
                            )
            

            for genoma in results:
                grafico_html = visualizza_geni_genoma(genoma)
                final_result.append({'genoma': genoma, 'grafico_html': grafico_html})

            return render(request, 'search_results.html', {'results': final_result})
    else:
        form = GenomeSearchForm()

    return render(request, 'search_form.html', {'form': form})
from .models import Profile, Genome, GeneProtein, AnnotationGenome, AnnotationProtein
from .forms import GenomeSearchForm,Upload_data, DownloadTextForm, ProteinAnnotate, SequenceProtein, GenomeAnnotate, SequenceGenome, UserRegistrationForm

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SearchIO

from io import StringIO

from django.shortcuts import render, redirect, redirect
from django.contrib.auth import authenticate, login, logout
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth.decorators import login_required
from django.contrib import messages

from django.utils import timezone
from django.db.models import Q
from django.http import StreamingHttpResponse

import re
import plotly.graph_objects as go

def home(request):

    status_genome = AnnotationGenome.objects.select_related("genome").filter(genome__sequence__isnull=False).order_by("annotation_time")[:10]
    status_proteins = AnnotationProtein.objects.select_related("geneprotein").order_by("annotation_time")[:64]

    print(status_genome)
    return render(request, 'home.html',{'status_genome':status_genome,'status_protein':status_proteins})

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
        form = UserRegistrationForm()

    return render(request, 'Authentication/register.html', {'form': form})


def login_view(request):
    if request.method == 'POST':
        form = AuthenticationForm(request, request.POST)
        if form.is_valid():
            username = form.cleaned_data.get('username')
            password = form.cleaned_data.get('password')
            user = authenticate(request, username=username, password=password, is_approved=True)

            if user is not None:
                if user.is_approved:
                    login(request, user)
                    next_url = request.POST.get('next', 'home')
                    if not next_url:
                        next_url = 'home'
                    print(next_url)
                    return redirect(next_url)
                else:
                    messages.error(request, 'Your account is not approved yet. Please wait for approval.')
            else:
                messages.error(request, 'Invalid username or password.')

    else:
        form = AuthenticationForm()

    return render(request, 'Authentication/login.html', {'form': form})


def logout_view(request):
        logout(request)
        return redirect('/')

def search_results(request):
    if request.method == 'POST':
        form = GenomeSearchForm(request.POST)
        if form.is_valid():
            sequence_query = form.cleaned_data['sequence']
            species = form.cleaned_data['species']
            chromosome = form.cleaned_data['chromosome']
            output_type = form.cleaned_data['output_type']
                   
            c = 0
            final_result = []

            if output_type == 'genome':

                results = AnnotationGenome.objects.select_related('genome').filter(
                                                                            Q(genome__is_validated=True) &
                                                                            Q(genome__chromosome__startswith=chromosome) &
                                                                            Q(species__contains=species))
            
            else:
                gene = form.cleaned_data['gene']
                transcript = form.cleaned_data['transcript']
                function = form.cleaned_data['function']
                 
                if validate(sequence_query,'dna'):
                   type_sequence = 'cds'        
                else:
                    type_sequence = 'pep'
                
                results = AnnotationProtein.objects.select_related('geneprotein__genome__annotationgenome').filter(
                                    Q(geneprotein__is_validated=True) &
                                    Q(geneprotein__type=type_sequence) &
                                    Q(gene__startswith=gene) &
                                    Q(description__contains=function) &
                                    Q(transcript__startswith=transcript) &
                                    Q(geneprotein__genome__chromosome__startswith=chromosome) &
                                    Q(geneprotein__genome__annotationgenome__species__contains=species)
                                    )       
            
            for DNAsequence in results:
                if output_type == 'genome':
                    my_dna = Seq(DNAsequence.genome.sequence)
                    if my_dna.count(sequence_query) > 0:
                        result_dic = {
                                        'type': 'genome',
                                        'chromosome': DNAsequence.genome.chromosome,
                                        'id': DNAsequence.genome.id,
                                    }
                else:
                    my_protein = Seq(DNAsequence.geneprotein.sequence)
                    if my_protein.count(sequence_query):
                        result_dic = {
                                        'type': type_sequence,
                                        'chromosome': DNAsequence.geneprotein.accession_number,
                                        'id': DNAsequence.geneprotein.id,
                                    }
                
                final_result.append(result_dic)
                print(final_result)
            
            return render(request, 'Search/search_results.html', {'results': final_result})
        else:
            return render(request, 'Search/search_form.html', {'form': form})
    else:
        return render(request, 'Search/search_form.html', {'form': GenomeSearchForm})

def validate(seq, alphabet='dna'):
    
    alphabets = {'dna': re.compile('^[acgtn]*$', re.I), 
             'protein': re.compile('^[acdefghiklmnpqrstvwy]*$', re.I)}


    if alphabets[alphabet].search(seq) is not None:
         return True
    else:
         return False

@login_required(login_url="/login")
def upload(request):
    if request.POST:
        form = Upload_data(request.POST, request.FILES)
        if form.is_valid():
            file = request.FILES['sequence']
            if 'annotated' not in request.POST:
                annotated = False
            else:
                annotated = request.POST['annotated']
            type = request.POST['output_type']
            file_content = file.read().decode('utf-8')
            file_name = file.name 
            file_stream = StringIO(file_content)
            username_id = request.user.id

            if type == 'genome':
                
                try:
                    genome = SeqIO.read(file_stream, 'fasta')
                except:
                    error_message = 'Genome expects exactly one fasta sequence'
                    return render(request, 'upload.html', {'form': form, 'error_message': error_message})

                sequence = genome.seq
                description = genome.description.split(':')
                chromosome = description[2]

                if not sequence:
                        error_message = 'The genome '+ chromosome + ' does not contain a sequence'
                        return render(request, 'upload.html', {'form': form, 'error_message': error_message})

                start = description[4]
                end = description[5]
                    
                genome = Genome.objects.get_or_create(chromosome=chromosome,
                                                    defaults={
                                                        'sequence':sequence,
                                                        'start':start,
                                                        'end':end,
                                                        'upload_time':timezone.now()
                                                            })
                                                
                # We populate the "artificial" genome created by the proteins
                if not genome[1] and genome[0].sequence is None:
                    genome[0].sequence=sequence
                    genome[0].start=start
                    genome[0].end=end
                    genome[0].save()
                
                elif genome[1]:
                    genome[0].save()
                
                else:
                    error_message = 'This genome is already in the database'
                    return render(request, 'upload.html', {'form': form, 'error_message': error_message})
                               
                if annotated:
                    species = re.split(r'[_.]', file_name)[:-1]
                    species = " ".join(species)
                    id_user = Profile.objects.get(id=username_id)
                    annotations = AnnotationGenome.objects.create(species=species,
                                                                    genome=genome[0],
                                                                    annotator=id_user,
                                                                    is_annotated=True,
                                                                    annotation_time = timezone.now())
                    annotations.save()
                return redirect(home)
            
            else :
                for record in SeqIO.parse(file_stream, 'fasta'):
                    a_n = record.id
                    sequence = record.seq

                    if not sequence:
                        error_message = 'The protein '+ a_n + ' does not contain a sequence'
                        return render(request, 'upload.html', {'form': form, 'error_message': error_message})

                    type = record.description.split()[1].lower()
                    description = re.split(r'[: ]', record.description)
                    start = description[5]
                    end = description[6]
                    chromosome = description[3]


                    genome, created = Genome.objects.get_or_create(chromosome=chromosome,
                                                    defaults={
                                                        'upload_time':timezone.now()
                                                            })

                    if created and genome.sequence is None:
                        genome.save()
                   
                    gene_protein, created_protein = GeneProtein.objects.get_or_create(accession_number=a_n,
                                                        defaults={
                                                            'sequence':sequence,
                                                            'type':type,
                                                            'start':start,
                                                            'end':end,
                                                            'genome':genome,
                                                            'upload_time':timezone.now()
                                                            })
                    
                    if created_protein:
                        gene_protein.save()
                    else:
                        error_message = 'The protein '+ a_n + ' is already in the database'
                        return render(request, 'upload.html', {'form': form, 'error_message': error_message})



                    if annotated and created_protein:
                        print(description)
                        gene=description[9]
                        transcript=a_n
                        gene_biotype=description[12]
                        gene_symbol=description[15]
                        transcript_biotype=description[14]

                        if type == "cds":
                            description_protein=' '.join(description[18:21])
                        else:
                            description_protein=' '.join(description[19:21])

                        id_user = Profile.objects.get(id=username_id)

                        annotation_protein = AnnotationProtein.objects.create(
                                                                        gene=gene,
                                                                        transcript=transcript,
                                                                        gene_biotype=gene_biotype, 
                                                                        transcript_biotype=transcript_biotype,
                                                                        gene_symbol=gene_symbol, 
                                                                        description=description_protein, 
                                                                        annotator=id_user,
                                                                        geneprotein=gene_protein,
                                                                        is_annotated=True,
                                                                        annotation_time = timezone.now()) 
                                                                    
                        annotation_protein.save()

                return redirect(home)
        else:
            return render(request,'upload.html',{'form':form})

    return render(request,'upload.html',{'form':Upload_data})


def visualisation_sequence(request,type,id):
    if type == "genome":
        genome = Genome.objects.get(id=id)
    else:
        genome = GeneProtein.objects.get(id=id,type=type)
    return render(request, 'Search/visualisation_sequence.html', {'genome_id': genome})


@login_required(login_url="/login")
def validator_view(request):
    unvalidated_genomes = Genome.objects.filter(is_validated=False, annotationgenome__isnull=True)
    unvalidated_proteins = GeneProtein.objects.filter(is_validated=False,annotationprotein__isnull=True)
    annotators = None


    if unvalidated_genomes != None or unvalidated_proteins != None:
        annotators = Profile.objects.filter(role='annotator')
    
        

    context = {
        'annotators': annotators,
        'unvalidated_genomes': unvalidated_genomes,
        'unvalidated_proteins': unvalidated_proteins
    }

    return render(request, 'Validator/validator_view.html', context)

@login_required(login_url="/login")
def validate_include_database(request):
    id_sequence = request.session.get('id_sequence')
    type_sequence = request.session.get('type')

    if request.method == 'POST':
        if type_sequence == "genome":
            form_annotate = GenomeAnnotate(request.POST)
            form_sequence = SequenceGenome(request.POST)
        else:
            form_annotate = ProteinAnnotate(request.POST)
            form_sequence = SequenceProtein(request.POST)

        if form_annotate.is_valid() and form_sequence.is_valid():
                if type_sequence == "genome":
                    validate_sequence = Genome.objects.get(pk=id_sequence)
                else:
                    validate_sequence = GeneProtein.objects.get(pk=id_sequence)
            
                validate_sequence.is_validated = True
                validate_sequence.save()

                return render(request, 'status_sequence_modified.html')

        else:
                return render(request,'Validator/validate_sequences.html',{'form': form_annotate , 'form2': form_sequence})

    return render(request,'Validator/validate_sequences.html',{'form': form_annotate , 'form2': form_sequence})
        

@login_required(login_url="/login")
def assigned_annotators(request):
    if request.method == 'POST':
   
        anno_id_protein = request.POST.getlist('annotation_id_protein')
        anno_id_genome = request.POST.getlist('annotation_id_genome')

        annotator_protein = request.POST.getlist('annotator_protein')
        annotator_genome = request.POST.getlist('annotator_genome')
        
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



@login_required(login_url="/login")
def annotator_view(request):
    annotator = request.user.id
    unannotated_genomes = AnnotationGenome.objects.filter(is_annotated=False, 
                                                          annotator=annotator,
                                                          genome__sequence__isnull=False)
    unannotated_proteins= AnnotationProtein.objects.filter(is_annotated=False, annotator=annotator)[:6]

    context = {
        'protein_annotations': unannotated_proteins,
        'genome_annotations': unannotated_genomes,
    }

    return render(request, 'Annotator/Annotator_dashboard.html', context)

@login_required(login_url="/login")
def sequence_view(request,type_of_sequence,id):
    request.session['type'] = type_of_sequence
    if type_of_sequence == "genome":
        print(id)
        genome_annotate_instance = AnnotationGenome.objects.get(pk=id)
        sequence_genome_instance = Genome.objects.get(pk=genome_annotate_instance.genome.id)
        
        request.session['id_annotate'] = id
        request.session['id_sequence'] = sequence_genome_instance.id

        annotate_form = GenomeAnnotate(instance=genome_annotate_instance)
        sequence_form = SequenceGenome(instance=sequence_genome_instance)
    else:

        protein_annotate_instance = AnnotationProtein.objects.get(pk=id)
        sequence_protein_instance = GeneProtein.objects.get(pk=protein_annotate_instance.geneprotein.id)
        
        request.session['id_annotate'] = id
        request.session['id_sequence'] = sequence_protein_instance.id

        annotate_form  = ProteinAnnotate(instance=protein_annotate_instance)
        sequence_form = SequenceProtein(instance=sequence_protein_instance)
    
    if "annotation" in request.path:
        return render(request,'Annotator/Annotate_sequences.html',{'form': annotate_form , 'form2': sequence_form})
    else:
        print('hey')
        return render(request,'Validator/validate_sequences.html',{'form': annotate_form , 'form2': sequence_form})

def save_annotation(request):
    id_annotate = request.session.get('id_annotate')
    id_sequence = request.session.get('id_sequence')
    type_sequence = request.session.get('type')

    if request.method == 'POST':

        if type_sequence == "genome":
            form_annotate = GenomeAnnotate(request.POST)
            form_sequence = SequenceGenome(request.POST)
        else:
            form_annotate = ProteinAnnotate(request.POST)
            form_sequence = SequenceProtein(request.POST)

        
        if form_annotate.is_valid() and form_sequence.is_valid():
                annotate = form_annotate.cleaned_data
                sequence = form_sequence.cleaned_data

                if type_sequence == "genome":
                    update_annotation = AnnotationGenome.objects.get(pk=id_annotate)
                    update_sequence = Genome.objects.get(pk=id_sequence)
                else:
                    update_annotation = AnnotationProtein.objects.get(pk=id_annotate)
                    update_sequence = GeneProtein.objects.get(pk=id_sequence)
                    

                for field_name, field_value in annotate.items():
                    setattr(update_annotation, field_name, field_value)

                update_annotation.is_annotated = True
                update_annotation.annotation_time = timezone.now()
                update_annotation.save()
        
                for field_name, field_value in sequence.items():
                    setattr(update_sequence, field_name, field_value)
                
                update_sequence.save()

                return render(request, 'status_sequence_modified.html')

        else:
                return render(request,'Annotator/Annotate_sequences.html',{'form': form_annotate , 'form2': form_sequence})

    return render(request,'Annotator/Annotate_sequences.html',{'form': form_annotate , 'form2': form_sequence})

def text_extraction(request):
    if request.method == 'POST':
        form = DownloadTextForm(request.POST)
        if form.is_valid():
            startposition = form.cleaned_data['start_position']
            endposition = form.cleaned_data['end_position']
            chromosome = form.cleaned_data['chromosome']
            species_query = form.cleaned_data['species']
            output_type = form.cleaned_data['output_type']

            gene = form.cleaned_data['gene']
            transcript = form.cleaned_data['transcript']
            function = form.cleaned_data['function']


            if output_type == 'genome':
                results = AnnotationGenome.objects.select_related('genome').filter(Q(genome__is_validated = True) &
                                                                                   (Q(species__contains=species_query) &
                                                                                   Q(genome__chromosome__startswith=chromosome)))
            else:
                results = AnnotationProtein.objects.select_related('geneprotein__genome__annotationgenome').filter(
                                                                            Q(geneprotein__is_validated=True) &
                                                                            Q(gene__startswith=gene)&
                                                                            Q(transcript__startswith=transcript) &
                                                                            Q(description__contains=function) &
                                                                            (Q(geneprotein__genome__chromosome__startswith=chromosome) &
                                                                            Q(geneprotein__genome__annotationgenome__isnull=False) & 
                                                                            Q(geneprotein__genome__annotationgenome__species__contains=species_query))
                                                                        )
        
            lines = []
            flag_None = False

            if startposition is None:
                startposition = 0

            if endposition is None:
                flag_None = True
            
            species = 0

            if not results.exists():
                lines.append(f'No results') 

            elif output_type == "genome":
                for type_sequence in results:
                    species += 1

                    my_sequence = Seq(type_sequence.genome.sequence)

                    if flag_None:
                        endposition = len(my_sequence)
                
                    print(startposition, endposition)

                    f = SeqFeature(FeatureLocation(startposition, endposition), type="CDS")
                    if f:
                        extraction = f.extract(my_sequence)
                        lines.append(f'>{type_sequence.species};{type_sequence.genome.chromosome};{startposition};{endposition}\n')
                        
                        for i in range(0, len(extraction), 70):
                            lines.append(f'{extraction[i:i+70]}\n')
            else:
                for type_sequence in results:

                    species += 1

                    my_sequence = Seq(type_sequence.geneprotein.sequence)
                    these_species = type_sequence.geneprotein.genome.annotationgenome.species
                    this_gene = type_sequence.gene
                    this_function = type_sequence.description
                    this_transcript = type_sequence.transcript

                    lines.append(f'>{these_species};{type_sequence.geneprotein.type};{this_gene};{this_transcript};{this_function}\n')
                            
                    for i in range(0, len(my_sequence), 70):
                                lines.append(f'{my_sequence[i:i+70]}\n')

                  
            response = StreamingHttpResponse((line for line in lines), content_type='text/plain')
            response['Content-Disposition'] = f'attachment; filename={species}_{output_type}_found.txt'
            return response
    else:
        form = DownloadTextForm()

    return render(request, 'text_extraction.html', {'form': form})




def visualizza_geni_genoma(genoma):
    # Recupera tutti i geni annotati sul genoma specificato
    geni_annotati = GeneProtein.objects.filter(annotated=False, genome=genoma)

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

    # Aggiungi tutte le tracce al layout
    layout = go.Layout(title='Geni sul genoma',
                       xaxis=dict(title='Posizione nel genoma'),
                       yaxis=dict(visible=False),
                       showlegend=True)

    # Crea la figura
    figura = go.Figure(data=tracce, layout=layout)

    # Converti la figura in HTML
    grafico_html = figura.to_html(full_html=False)

    return grafico_html

def visualizzazione(request):
    # Recupera un'istanza di GeneProtein (esempio)
    geneprotein_instance = GeneProtein.objects.first()

    # Ora otteniamo l'istanza di Genome associata a questa istanza di GeneProtein
    genome_instance = geneprotein_instance.genome

    # Ora possiamo chiamare la funzione visualizza_geni_genoma con l'istanza di Genome
    grafico_html = visualizza_geni_genoma(genome_instance)

    return render(request, 'visualization.html', {'grafico_html': grafico_html})

def blast_search(request):
    # Replace 'your_sequence' with your actual DNA, RNA, or protein sequence
    ## récupérer la séquence mise dans query 
    # if request.method == 'POST':
    #     form = GenomeSearchForm(request.POST)
    #     if form.is_valid():
    #         if form.cleaned_data['database'] == 'NCBI_Blast':
    #             sequence_query = form.cleaned_data['sequence']
    # else:
    #     form = GenomeSearchForm()
        
    sequence = "AAAAGTGTACGGATTCTGGAAGCTGAATGCTGTGCAGATCATATCCATATGCTTGTGGAG"
    # Specify the BLAST program (e.g., 'blastn' for nucleotides, 'blastp' for proteins)
    blast_program = "blastn" ## faire des conditions selon ce que l'utilisateur a mis comme seq (peptide ou ADN)
    try:
        result_handle = NCBIWWW.qblast(program=blast_program, database="nt", sequence=sequence, descriptions=50, hitlist_size=25)
        blast_results = SearchIO.read(result_handle, "blast-xml")
        context = {'blast_results': blast_results}
        return render(request, 'blast_results.html', context)
    except Exception as e:
        return "An error occurred"
from .models import Profile, Genome, GeneProtein, AnnotationGenome, AnnotationProtein, AnnotationStatus
from .forms import GenomeSearchForm,Upload_data, DownloadTextForm, ProteinAnnotate, SequenceProtein, GenomeAnnotate, SequenceGenome, UserRegistrationForm, CommentForm
from django.contrib.contenttypes.models import ContentType
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

from django.http import HttpResponseForbidden
from django.urls import reverse
from django.http import HttpResponseRedirect

from django.utils import timezone
from django.db.models import Q
from django.http import StreamingHttpResponse

import re
import plotly.graph_objects as go
from functools import wraps

from django.shortcuts import get_object_or_404

from django.core.paginator import Paginator


def home(request):
    return render(request, 'home.html')

def user_registration_login_status(view_func):
    @wraps(view_func)
    def _wrapped_view(request, *args, **kwargs):
        if ("register" in request.path or "login" in request.path) and request.user.is_authenticated:
            return HttpResponseRedirect(reverse('home'))
        else:
            return view_func(request, *args, **kwargs)
    
    return _wrapped_view

def user_is_validator(view_func):
   
    @wraps(view_func)
    def _wrapped_view(request, *args, **kwargs):

        if (request.user.role == 'validator' or request.user.is_superuser) and request.user.is_authenticated:
            return view_func(request, *args, **kwargs)
        else:
            return HttpResponseForbidden(f"Access Forbidden - {request.user.username} is not validator or admin")
    return _wrapped_view

def user_is_annotator(view_func):
   
    @wraps(view_func)
    def _wrapped_view(request, *args, **kwargs):

        if (request.user.role == 'annotator' or request.user.is_superuser) and request.user.is_authenticated:
            return view_func(request, *args, **kwargs)
        else:
            return HttpResponseForbidden(f"Access Forbidden - {request.user.username} is not annotator or admin")
    return _wrapped_view

def is_not_user(view_func):
   
    @wraps(view_func)
    def _wrapped_view(request, *args, **kwargs):

        if request.user.role != 'user'  and request.user.is_authenticated:
            return view_func(request, *args, **kwargs)
        else:
            return HttpResponseForbidden(f"Access Forbidden - {request.user.username} is not an annotator, validator or admin")
    return _wrapped_view


@user_registration_login_status
def register_view(request):
    if request.method == 'POST':
        form = UserRegistrationForm(request.POST)
        if form.is_valid():
            user = form.save(commit=False)
            user.username = user.username.lower()
            user.save()
            return render(request,'home.html')
    else:
        form = UserRegistrationForm()

    return render(request, 'Authentication/register.html', {'form': form})


@user_registration_login_status
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
    if request.method == 'GET' and request.GET:
        form = GenomeSearchForm(request.GET)
        if form.is_valid():
            sequence_query = form.cleaned_data['sequence']
            species = form.cleaned_data['species']
            chromosome = form.cleaned_data['chromosome']
            output_type = form.cleaned_data['output_type']
            database = form.cleaned_data['database']

            final_result = []

            if database == "BactaHub":
                sequence_type = "unknown"
                if output_type == 'genome':

                    results = AnnotationGenome.objects.select_related('genome').filter(
                                                                                Q(genome__is_validated=True) &
                                                                                Q(genome__chromosome__startswith=chromosome) &
                                                                                Q(species__contains=species)).only('genome__chromosome','species')
                
                else:
                    gene = form.cleaned_data['gene']
                    transcript = form.cleaned_data['transcript']
                    function = form.cleaned_data['function']
                 
                    if validate(sequence_query,'dna'):
                        sequence_type = 'cds'        
                    else:
                        sequence_type = 'pep'
                
                    results = AnnotationProtein.objects.select_related('geneprotein__genome__annotationgenome').filter(
                                        Q(geneprotein__is_validated=True) &
                                        Q(geneprotein__type=sequence_type) &
                                        Q(gene__startswith=gene) &
                                        Q(description__contains=function) &
                                        Q(transcript__startswith=transcript) &
                                        Q(geneprotein__genome__chromosome__startswith=chromosome) &
                                        Q(geneprotein__genome__annotationgenome__species__contains=species)
                                        ).only('geneprotein__type','geneprotein__accession_number','geneprotein__sequence','geneprotein__genome__id','geneprotein__genome__annotationgenome__species')       
                
                special_query = False

                if re.search('%',sequence_query):
                    special_query = True
                    sequence_pattern = sequence_query.split("%")
                    if sequence_type == 'cds' or output_type == 'genome':
                        pattern = re.compile(sequence_pattern[1] + "([ATGCN])*?" + sequence_pattern[2])
                    else:
                        pattern = re.compile(sequence_pattern[1] + "([ACDEFGHIKLMNPQRSTVWY])*?" + sequence_pattern[2])


                if special_query:
                    print(pattern)
                    for DNAsequence in results:
                        if output_type == 'genome': 
                            if pattern.search(DNAsequence.genome.sequence):
                                result_dic = {
                                                'type': 'genome',
                                                'species': DNAsequence.species,
                                                'chromosome': DNAsequence.genome.chromosome,
                                                'id': DNAsequence.genome.id,
                                            }
                                final_result.append(result_dic)
                        else:
                            if pattern.search(DNAsequence.geneprotein.sequence):
                                
                                result_dic = {
                                                'type': sequence_type,
                                                'species': DNAsequence.geneprotein.genome.annotationgenome.species,
                                                'chromosome': DNAsequence.geneprotein.accession_number,
                                                'id': DNAsequence.geneprotein.id,
                                            }
                        
                                final_result.append(result_dic)

                else:
                    for DNAsequence in results:
                        if output_type == 'genome':
                            my_dna = Seq(DNAsequence.genome.sequence)
                            if my_dna.count(sequence_query) > 0:
                                result_dic = {
                                                'type': 'genome',
                                                'species': DNAsequence.species,
                                                'chromosome': DNAsequence.genome.chromosome,
                                                'id': DNAsequence.genome.id,
                                            }
                                final_result.append(result_dic)
                        else:
                            my_protein = Seq(DNAsequence.geneprotein.sequence)
                            if my_protein.count(sequence_query):
                                result_dic = {
                                                'type': sequence_type,
                                                'species': DNAsequence.geneprotein.genome.annotationgenome.species,
                                                'chromosome': DNAsequence.geneprotein.accession_number,
                                                'id': DNAsequence.geneprotein.id,
                                            }
                        
                                final_result.append(result_dic)

                p = Paginator(final_result,15)
                page = request.GET.get('page')
                sequences = p.get_page(page)

                return render(request, 'Search/search_results.html', {'results': sequences, 'form':form})
            
            else :
                if validate(sequence_query) == True: ## ADN
                    sequence = sequence_query
                    database = "nt"
                    program = "blastn"
                    return blast_search(request, sequence = sequence, database = database, program = program)

                elif validate(sequence_query) == False : ## protein
                    sequence = sequence_query
                    database = "uniprotkb"
                    program = "blastp"
                    return blast_search(request, sequence = sequence, database = database, program = program)

                else :
                    return render(request, 'Search/search_form.html', {'form': form})
                
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
                   
                    gene_protein, created_protein = GeneProtein.objects.get_or_create(accession_number=a_n,type=type,
                                                        defaults={
                                                            'sequence':sequence,
                                                            'start':start,
                                                            'end':end,
                                                            'genome':genome,
                                                            'upload_time':timezone.now()
                                                            })
                    
                    if not created_protein:
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
@user_is_validator
def validator_view(request):


    created_genomes = Genome.objects.filter(Q(is_validated=False) &
                                            Q(annotationgenome__isnull=True)).defer('sequence')

    created_proteins = GeneProtein.objects.filter(Q(is_validated=False) &
                                                        Q(annotationprotein__isnull=True)).defer('sequence')
    
    unvalidated_genomes = AnnotationGenome.objects.select_related('genome').filter(genome__is_validated=False, is_annotated=True)
                                           
    unvalidated_proteins = AnnotationProtein.objects.select_related('geneprotein').filter(geneprotein__is_validated=False, is_annotated=True)
    annotators = None


    if created_genomes != None or created_proteins != None:
        annotators = Profile.objects.filter(role='annotator')
    
    context = {
        'annotators': annotators,
        'created_genomes': created_genomes,
        'created_proteins': created_proteins,
        'unvalidated_genomes': unvalidated_genomes,
        'unvalidated_proteins': unvalidated_proteins
    }

    return render(request, 'Validator/validator_view.html', context)

@login_required(login_url="/login")
@user_is_validator
def validate_include_database(request, sequence_type=None, id_get_sequence=None):
    if request.method == 'POST':
            id_sequence = request.session.get('id_sequence')
            sequence_type = request.session.get('type')
            id_get_sequence = request.session.get('id_annotate')
            print(request.POST)

            if sequence_type == "genome":
                    form_annotate = GenomeAnnotate(request.POST, current_user=request.user)
                    form_sequence = SequenceGenome(request.POST, current_user=request.user)
            else:
                    form_annotate = ProteinAnnotate(request.POST, current_user=request.user)
                    form_sequence = SequenceProtein(request.POST, current_user=request.user)

            annotation_status = CommentForm(request.POST)

            if form_annotate.is_valid() and form_sequence.is_valid() and annotation_status.is_valid():  
                
                status_of_annotation = annotation_status.cleaned_data['status']
                comment_for_annotator = annotation_status.cleaned_data['comment']
                
                if sequence_type == "genome":
                        validate_sequence = Genome.objects.get(pk=id_sequence)
                        sequence = form_sequence.cleaned_data['chromosome']
                        type_object = ContentType.objects.get_for_model(AnnotationGenome)
                        annotation = AnnotationGenome.objects.get(pk=id_get_sequence)
                        
                else:   
                        type_object = ContentType.objects.get_for_model(AnnotationProtein)
                        annotation = AnnotationProtein.objects.get(pk=id_get_sequence)
                        validate_sequence = GeneProtein.objects.get(pk=id_sequence)
                        sequence = form_sequence.cleaned_data['accession_number']

                if status_of_annotation == "validated":       
                    validate_sequence.is_validated = True
                    validate_sequence.save()
                
                else:
                    annotation.is_annotated = False
                    annotation.save()

                if comment_for_annotator:
                    validator = Profile.objects.get(pk = request.user.id)

                    AnnotationStatus.objects.create(content_type=type_object,
                                                            object_id = id_get_sequence,
                                                            validator = validator,
                                                            comment = comment_for_annotator,
                                                            status = status_of_annotation)

                messages.success(request, f'Sequence "{sequence}" has been {status_of_annotation}')
                return HttpResponseRedirect(reverse('validator_view'))

            else:
                return render(request,'Validator/validate_sequences.html',{'form': form_annotate , 'form2': form_sequence, 'comment': annotation_status})
    if request.method == 'GET':
        print(sequence_type)
        if sequence_type == "genome":
            validate_sequence = get_object_or_404(AnnotationGenome,pk=id_get_sequence)
            validate_sequence.genome.is_validated = True
            validate_sequence.genome.save()
            sequence = validate_sequence.genome.chromosome
            messages.success(request, f'Sequence "{sequence}" has been included in the database')
                        
        else:   
            validate_sequence = get_object_or_404(AnnotationProtein,pk=id_get_sequence)
            validate_sequence.geneprotein.is_validated = True
            validate_sequence.geneprotein.save()
            sequence = validate_sequence.geneprotein.accession_number
            messages.success(request, f'Sequence "{sequence}" has been included in the database')
        
        return HttpResponseRedirect(reverse('validator_view'))
    
    return render(request,'Validator/validate_sequences.html',{'form': form_annotate , 'form2': form_sequence, 'comment': annotation_status})
        

@login_required(login_url="/login")
@user_is_validator
def assigned_annotators(request):
    if request.method == 'POST':

        print(request.POST)
   
        anno_id_protein = request.POST.getlist('annotation_id_protein')
        anno_id_genome = request.POST.getlist('annotation_id_genome')


        annotator_protein = request.POST.getlist('annotator_protein')
        annotator_genome = request.POST.getlist('annotator_genome')
        
        if len(anno_id_protein) != 0:
            for i in range(len(anno_id_protein)):
                new_annotation_protein = AnnotationProtein.objects.create(
                    annotator_id=annotator_protein[i],
                    geneprotein_id=anno_id_protein[i]
                )
                new_annotation_protein.save()
        if len(anno_id_genome) != 0:
            for i in range(len(anno_id_genome)):
                new_annotation_genome = AnnotationGenome.objects.create(
                    annotator_id=annotator_genome[i],
                    genome_id=anno_id_genome[i]
                )
                new_annotation_genome.save()

        return render(request, 'Validator/validate_annotators.html')

    return redirect('home')

@login_required(login_url="/login")
@user_is_validator
def delete_sequence_from_database(request,type_of_sequence,id):
    request.session['type'] = type_of_sequence
    if type_of_sequence == "genome":
        genome_annotate_instance = get_object_or_404(AnnotationGenome, pk=id)
        sequence = genome_annotate_instance.genome.chromosome
        genome_annotate_instance.genome.delete()
        
    else:
        protein_annotate_instance = get_object_or_404(AnnotationProtein, pk=id)
        sequence = protein_annotate_instance.geneprotein.accession_number
        protein_annotate_instance.geneprotein.delete()
    
    messages.info(request, f'Sequence "{sequence}" has been deleted successfully.')
    return HttpResponseRedirect(reverse('validator_view'))


@login_required(login_url="/login")
@user_is_annotator
def annotator_view(request):
    annotator = request.user.id
    unannotated_genomes = AnnotationGenome.objects.filter(is_annotated=False, 
                                                          annotator=annotator,
                                                          genome__sequence__isnull=False)
    unannotated_proteins= AnnotationProtein.objects.filter(is_annotated=False, annotator=annotator)[:6]

    status_genome = AnnotationGenome.objects.select_related('genome').filter(annotator_id = request.user.id)
    status_protein = AnnotationProtein.objects.select_related('geneprotein').filter(annotator_id = request.user.id)
                                                           
    context = {
        'protein_annotations': unannotated_proteins,
        'genome_annotations': unannotated_genomes,
        'status_genome': status_genome,
        'status_protein': status_protein
    }

    return render(request, 'Annotator/Annotator_dashboard.html', context)

@login_required(login_url="/login")
@is_not_user
def sequence_view(request,type_of_sequence,id):
    request.session['type'] = type_of_sequence
    print('hey')
    if type_of_sequence == "genome":
        genome_annotate_instance = AnnotationGenome.objects.get(pk=id)
        sequence_genome_instance = Genome.objects.get(pk=genome_annotate_instance.genome.id)
        
        request.session['id_annotate'] = id
        request.session['id_sequence'] = sequence_genome_instance.id

        annotate_form = GenomeAnnotate(instance=genome_annotate_instance, current_user=request.user)
        sequence_form = SequenceGenome(instance=sequence_genome_instance, current_user=request.user)
    else:

        protein_annotate_instance = AnnotationProtein.objects.get(pk=id)
        sequence_protein_instance = GeneProtein.objects.get(pk=protein_annotate_instance.geneprotein.id)
        
        request.session['id_annotate'] = id
        request.session['id_sequence'] = sequence_protein_instance.id

        annotate_form  = ProteinAnnotate(instance=protein_annotate_instance, current_user=request.user)
        sequence_form = SequenceProtein(instance=sequence_protein_instance, current_user=request.user)
    
    if "annotation" in request.path:
        return render(request,'Annotator/Annotate_sequences.html',{'form': annotate_form , 'form2': sequence_form})
    else:
        comment_form = CommentForm()
        return render(request,'Validator/validate_sequences.html',{'form': annotate_form , 'form2': sequence_form, 'comment': comment_form})

@login_required(login_url="/login")
@user_is_annotator
def save_annotation(request):
    id_annotate = request.session.get('id_annotate')
    id_sequence = request.session.get('id_sequence')
    sequence_type = request.session.get('type')

    if request.method == 'POST':

        if sequence_type == "genome":
            form_annotate = GenomeAnnotate(request.POST, current_user=request.user)
            form_sequence = SequenceGenome(request.POST, current_user=request.user)
        else:
            form_annotate = ProteinAnnotate(request.POST, current_user=request.user)
            form_sequence = SequenceProtein(request.POST, current_user=request.user)

        print(form_sequence)

        if form_annotate.is_valid() and form_sequence.is_valid():
                
                annotate = form_annotate.cleaned_data
                sequence = form_sequence.cleaned_data

                if sequence_type == "genome":
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

@login_required(login_url="/login")
def text_extraction(request, sequence_type = None, id_sequence = None):
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
                for sequence_type in results:
                    species += 1

                    my_sequence = Seq(sequence_type.genome.sequence)

                    if flag_None:
                        endposition = len(my_sequence)
                
                    print(startposition, endposition)

                    f = SeqFeature(FeatureLocation(startposition, endposition), type="CDS")
                    if f:
                        extraction = f.extract(my_sequence)
                        lines.append(f'>{sequence_type.species};{sequence_type.genome.chromosome};{startposition};{endposition}\n')
                        
                        for i in range(0, len(extraction), 70):
                            lines.append(f'{extraction[i:i+70]}\n')
            else:
                for sequence_type in results:

                    species += 1

                    my_sequence = Seq(sequence_type.geneprotein.sequence)
                    these_species = sequence_type.geneprotein.genome.annotationgenome.species
                    this_gene = sequence_type.gene
                    this_function = sequence_type.description
                    this_transcript = sequence_type.transcript

                    lines.append(f'>{these_species};{sequence_type.geneprotein.type};{this_gene};{this_transcript};{this_function}\n')
                            
                    for i in range(0, len(my_sequence), 70):
                                lines.append(f'{my_sequence[i:i+70]}\n')

                  
            response = StreamingHttpResponse((line for line in lines), content_type='text/plain')
            response['Content-Disposition'] = f'attachment; filename={species}_{output_type}_found.txt'
            return response
    
    elif request.method == 'GET' and "download" in request.get_full_path():
               
                output_type = sequence_type
    
                if output_type == 'genome':
                    results = AnnotationGenome.objects.select_related('genome').filter(Q(genome__is_validated = True) &
                                                                                    (Q(genome__pk = id_sequence)))
                else:
                    results = AnnotationProtein.objects.select_related('geneprotein__genome__annotationgenome').filter(
                                                                                Q(geneprotein__is_validated=True) &
                                                                                Q(geneprotein__pk = id_sequence))
            
                lines = []
            
                if not results.exists():
                    lines.append(f'No results') 

                elif output_type == "genome":
                    for sequence_type in results:
                            txt_title = sequence_type.genome.chromosome
                            extraction = sequence_type.genome.sequence
                            species = sequence_type.species
                            lines.append(f'>{species};{sequence_type.genome.chromosome};{sequence_type.genome.start};{sequence_type.genome.end}\n')
                            
                            for i in range(0, len(extraction), 70):
                                lines.append(f'{extraction[i:i+70]}\n')
                else:
                    for sequence_type in results:

                        txt_title = sequence_type.gene
                        my_sequence = sequence_type.geneprotein.sequence
                        species = sequence_type.geneprotein.genome.annotationgenome.species
                        this_gene = sequence_type.gene
                        this_function = sequence_type.description
                        this_transcript = sequence_type.transcript

                        lines.append(f'>{species};{sequence_type.geneprotein.accession_number};{sequence_type.geneprotein.type};{this_gene};{this_transcript};{this_function}\n')
                                
                        for i in range(0, len(my_sequence), 70):
                                    lines.append(f'{my_sequence[i:i+70]}\n')

                    
                response = StreamingHttpResponse((line for line in lines), content_type='text/plain')
                response['Content-Disposition'] = f'attachment; filename={species}_{output_type}_{txt_title}.txt'
                return response
            
    else:
        print(request.GET)
        form = DownloadTextForm()

    return render(request, 'text_extraction.html', {'form': form})

@login_required(login_url="/login")
def visualize_genome_genes(request, sequence_type, selected_genome_id):
 
    if sequence_type == "genome":
        annotated_genes = AnnotationProtein.objects.select_related('geneprotein').filter(geneprotein__is_validated = True, 
                                                                                         geneprotein__genome__id = selected_genome_id)[:1000]
    else:
        unknown_gene = AnnotationProtein.objects.select_related('geneprotein').get(geneprotein_id = selected_genome_id)
        annotated_genes = AnnotationProtein.objects.select_related('geneprotein').filter(Q(geneprotein__is_validated = True), 
                                                     Q(geneprotein__type = sequence_type), 
                                                     Q(geneprotein__genome__id = unknown_gene.geneprotein.genome.id),
                                                     Q(geneprotein__genome__id = unknown_gene.geneprotein.genome.id),
                                                     ~ Q(geneprotein__id = unknown_gene.id),
                                                     Q(geneprotein__start__gte = unknown_gene.geneprotein.start - 50000) &
                                                     Q(geneprotein__end__lte = unknown_gene.geneprotein.end + 50000))[:1000]
        annotated_genes = list(annotated_genes)  
        annotated_genes.append(unknown_gene)

    color_palette = ['rgb(31, 119, 180)', 'rgb(255, 127, 14)', 'rgb(44, 160, 44)', 'rgb(214, 39, 40)', 'rgb(148, 103, 189)',
                     'rgb(140, 86, 75)', 'rgb(227, 119, 194)', 'rgb(127, 127, 127)', 'rgb(188, 189, 34)', 'rgb(23, 190, 207)']

    traces = []
    for i, gene in enumerate(annotated_genes):
        if sequence_type != 'genome' and gene.id == unknown_gene.id:
            color = 'rgb(0, 0, 0)'  # Set a specific lor for the unknown protein
            line_width = 50
        else:
            color = color_palette[i % len(color_palette)]
            line_width = 15

        gene_trace = go.Scatter(
            x=[gene.geneprotein.start, gene.geneprotein.end],
            y=[0, 0],
            mode='lines',
            line=dict(color=color, width=line_width),
            hoverinfo='text',
            hovertext=f'Gene: {gene.gene} | Chromosome: {gene.geneprotein.genome.chromosome} | Start: {gene.geneprotein.start} | End: {gene.geneprotein.end}\
                        <br>Description: {gene.description}',
            name=gene.geneprotein.accession_number
        )
        traces.append(gene_trace)


    layout = go.Layout(title='Genes on the genome',
                       xaxis=dict(title='Position on the genome'),
                       yaxis=dict(visible=False),
                       showlegend=True,
                       width=1500,  
                       height=400,
                       paper_bgcolor='#fbfbf9')
    figure = go.Figure(data=traces, layout=layout)

 
    html_graph = figure.to_html(full_html=False)

    return render(request, 'visualization.html', {'html_graph': html_graph})



def blast_search(request,sequence,program,database):
        
    ##sequence = "AAAAGTGTACGGATTCTGGAAGCTGAATGCTGTGCAGATCATATCCATATGCTTGTGGAG"
    ##print("error 1")
    # Specify the BLAST program (e.g., 'blastn' for nucleotides, 'blastp' for proteins)
    ##blast_program = "blastn" ## faire des conditions selon ce que l'utilisateur a mis comme seq (peptide ou ADN)
    ##print("error 2")

    result_handle = NCBIWWW.qblast('blastn', 'nt', 'MT747438')

    try:
        result_handle = NCBIWWW.qblast(program = program, database= database, sequence=sequence, descriptions=50, hitlist_size=25)
        blast_results = SearchIO.read(result_handle, "blast-xml")
        context = {'blast_results': blast_results}
        print(blast_results)
        print("error 3")
        return render(request, 'Search/blast_results.html', context)
    except Exception as e:
        print("error 4")
        return "An error occurred"

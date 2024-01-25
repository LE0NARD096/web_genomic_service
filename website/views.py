from django.shortcuts import render,redirect
from .models import Profile, Genome, GeneProtein
from .forms import GenomeSearchForm,Upload_data, DownloadTextForm
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
            output_type = form.cleaned_data['output_type']

           
            c = 0
            final_result = []

            if output_type == 'genome':
                results = Genome.objects.all()

            else:
                 if validate(sequence_query,'dna'):
                    results = GeneProtein.objects.filter(type='cds')
                 else:
                     results = GeneProtein.objects.filter(type='pep')
        
            

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
                    genome.save()
                else:
                    genome = Genome.objects.create(sequence=sequence,
                                               chromosome=chromosome,
                                               start=start,
                                               end=end)
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
                    genome = Genome.objects.create(chromosome=chromosome)
                    genome.save()
                                           
                   gene_protein = GeneProtein.objects.create(accession_number=a_n,
                                                        sequence=sequence,
                                                        type=type,
                                                        start=start,
                                                        end=end,
                                                        genome=genome)
                   gene_protein.save()

            return redirect(home)
        else:
            return render(request,'upload.html',{'form':form})

    return render(request,'upload.html',{'form':Upload_data})


def visualisation_sequence(request,type,id):
    print('iii')
    if type == "genome":
        genome = Genome.objects.get(id=id)
    else:
        genome = GeneProtein.objects.get(id=id,type=type)
    return render(request, 'visualisation_sequence.html', {'genome_id': genome})

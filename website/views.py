from django.shortcuts import render,redirect
from .models import Profile, Genome
from .forms import GenomeSearchForm,Upload_data, DownloadTextForm
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO
from django.contrib.auth import authenticate, login, logout
from django.contrib import messages
from django.http import HttpResponse
from Bio.SeqFeature import SeqFeature, FeatureLocation, SimpleLocation
from django.urls import reverse

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
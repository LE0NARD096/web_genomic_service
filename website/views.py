from django.shortcuts import render,redirect
from .models import Profile, Genome
from .forms import GenomeSearchForm,Upload_data
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO
from django.contrib.auth import authenticate, login, logout
from django.contrib import messages
from django.http import HttpResponse

from django.urls import reverse

from django.contrib.auth.forms import AuthenticationForm, UserCreationForm
from django.contrib.auth.decorators import login_required
from django.shortcuts import redirect
from .forms import UserRegistrationForm

def home(request):
    # Add the following lines to get the URLs for register and login pages
    register_url = reverse('register')  # Replace 'register' with the actual name of your register URL
    login_url = reverse('login')    # Replace 'login_url' with the actual name of your login URL

    return render(request, 'home.html', {'register': register_url, 'login': login_url})


@login_required
def logout_view(request):
    logout(request)
    return redirect('home')  # Redirect to home after logout

def register_view(request):
    user = request.user
    print(user)
    #if user.is_authenticated:
        #return HttpResponse(f'You are already authenticated as {user.email}')
    if request.method == 'POST':
        form = UserRegistrationForm(request.POST)
        print(form)
        if form.is_valid():
            user = form.save(commit=False)
            user.username = user.username.lower()
            user.save()
            print('iii')
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
            user = authenticate(request, username=username, password=password)
            if user is not None:
                login(request, user)
                return redirect('upload_file') 
    else:
        form = AuthenticationForm()
    
    return render(request, 'Authentication/login.html', {'form': form})

def logout_view(request):
    logout(request)
    messages.success(request, f'You have been logged out.')
    return redirect('home')


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
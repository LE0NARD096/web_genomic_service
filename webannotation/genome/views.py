from django.shortcuts import render, get_object_or_404, redirect
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.decorators import login_required
from .models import Annotation, Sequence
from .models import User
from django.contrib.auth import authenticate, login
from django.contrib.auth.forms import AuthenticationForm
from .forms import CustomUserCreationForm, AnnotationForm, RoleUpdateForm


def home(request):
    return render(request, 'home.html')

def login_view(request):
    if request.method == 'POST':
        form = AuthenticationForm(request, request.POST)
        if form.is_valid():
            username = form.cleaned_data.get('username')
            password = form.cleaned_data.get('password')
            user = authenticate(request, username=username, password=password)
            if user is not None:
                login(request, user)
                return redirect('home')  # Redirect to home page after login
            else:
                # Handle invalid login
                pass
    else:
        form = AuthenticationForm()
    return render(request, 'login.html', {'form': form})


def user_list(request):
    if not request.user.is_superuser:
        return redirect('some_other_view')  # Rediriger si l'utilisateur n'est pas administrateur

    users = User.objects.all()
    return render(request, 'user_list.html', {'users': users})

# User Registration
def register(request):
    if request.method == 'POST':
        form = CustomUserCreationForm(request.POST)
        if form.is_valid():
            form.save()
            return redirect('login_url')  # Redirect to login page after registration
    else:
        form = CustomUserCreationForm()
    return render(request, 'register.html', {'form': form})


# User Profile
@login_required
def profile(request):
    return render(request, 'profile.html')

# Annotation Creation
def create_annotation(request):
    if request.method == 'POST':
        form = AnnotationForm(request.POST)
        if form.is_valid():
            form.save()
            return redirect('some_url')
    else:
        form = AnnotationForm()
    return render(request, 'create_annotation.html', {'form': form})

# Annotation Editing
def edit_annotation(request, annotation_id):
    annotation = Annotation.objects.get(id=annotation_id)
    if request.method == 'POST':
        form = AnnotationForm(request.POST, instance=annotation)
        if form.is_valid():
            form.save()
            return redirect('some_url')
    else:
        form = AnnotationForm(instance=annotation)
    return render(request, 'edit_annotation.html', {'form': form, 'annotation': annotation})

# Validate Annotation
def validate_annotation(request, annotation_id):
    # Implementation
    pass

# List Sequences
def list_sequences(request):
    sequences = Sequence.objects.all()
    return render(request, 'list_sequences.html', {'sequences': sequences})

# View Annotation
def view_annotation(request, annotation_id):
    annotation = Annotation.objects.get(id=annotation_id)
    return render(request, 'view_annotation.html', {'annotation': annotation})


def update_user_role(request, user_id):
    user = get_object_or_404(User, pk=user_id)
    if request.method == 'POST':
        form = RoleUpdateForm(request.POST, instance=user)  # Updated to RoleUpdateForm
        if form.is_valid():
            form.save()
            return redirect('some_view_name')
    else:
        form = RoleUpdateForm(instance=user)  # Updated to RoleUpdateForm

    return render(request, 'update_user_role.html', {'form': form, 'user': user})

from django import forms
from django.contrib.auth.forms import UserCreationForm
from .models import CustomUser
from .models import Annotation
from .models import Role, User

# Custom User Registration Form
class CustomUserCreationForm(UserCreationForm):
    first_name = forms.CharField(required=True)
    last_name = forms.CharField(required=True)
    phone_number = forms.CharField(required=True)
    email = forms.EmailField(required=True)
    ROLE_CHOICES = [
    ('user', 'User'),
    ('annotator', 'Annotator'),
    ('validator', 'Validator'),
]
    role = forms.ChoiceField(choices=ROLE_CHOICES)

    class Meta:
        model = CustomUser
        fields = ('username', 'first_name', 'last_name', 'phone_number', 'email', 'password1', 'password2', 'role')

    def save(self, commit=True):
        user = super().save(commit=False)
        # ... other field assignments ...
        role_name = self.cleaned_data['role']
        # Fetch the Role instance
        role_instance = Role.objects.get(name=role_name)
        user.role = role_instance
        if commit:
            user.save()
        return user





# Annotation Form
class AnnotationForm(forms.ModelForm):
    class Meta:
        model = Annotation
        fields = ['text', 'sequence', 'isValidated']

class RoleUpdateForm(forms.ModelForm):
    role = forms.ModelChoiceField(queryset=Role.objects.all(), required=True)

    class Meta:
        model = User
        fields = ['role']

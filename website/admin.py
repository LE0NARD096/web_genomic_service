from django.contrib import admin
from .models import Profile, Genome

@admin.register(Profile)
class ProfileAdmin(admin.ModelAdmin):
    fields = (('username','is_approved'),'email','role')
    list_display = ('username','role','is_approved')
    ordering = ('-is_approved',)

admin.site.register(Genome)


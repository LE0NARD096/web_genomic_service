"""
URL configuration for webannotation project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/5.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.urls import path, include
from website import views
from django.contrib import admin
from website.views import  blast_search
from django.conf import settings
from django.conf.urls.static import static


urlpatterns = [
    path('', views.home, name='home'),  # Add this line for the home page
    path("admin/", admin.site.urls),
    path('', include('website.urls')),
    path('blast_results/', blast_search, name='blast_results'),
    path('tinymce/', include('tinymce.urls')),
    # ... other url patterns ...
] + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)


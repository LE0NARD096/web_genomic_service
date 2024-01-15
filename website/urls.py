from django.urls import path
from . import views
from .views import search_results

urlpatterns = [
    path('',views.home,name = "home"),
    path('search/', search_results, name='search_results'),
]

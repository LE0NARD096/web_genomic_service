from django.http import HttpResponseRedirect
from django.urls import resolve, reverse

class RoleMiddleware:
    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        # Get the name of the view being requested
        current_view = resolve(request.path_info).url_name

        # If the request is for the login or register page, do not redirect
        if current_view in ['login_url', 'register']:
            return self.get_response(request)

        # Redirect unauthenticated users to the login page
        if not request.user.is_authenticated:
            return HttpResponseRedirect(reverse('login_url'))

        # Role-based redirection logic for authenticated users
        if request.user.is_authenticated:
            user_role = request.user.role.name  # Assuming 'role' is a field in your User model
            # Add your logic here for role-based redirection
            # Example: if user_role == 'lecteur': ...

        return self.get_response(request)

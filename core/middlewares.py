from django.http import HttpResponseForbidden
from django.conf import settings


class AdminLoginMiddleware(object):
    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request, *args, **kwargs):
        response = self.get_response(request)
        return response

    def process_view(self, request, view_func, *view_args, **view_kargs):
        # Get the view name as a string
        view_name = view_func.__module__

        # If the view name is in our exclusion list, exit early
        inclusion_set = getattr(settings, 'INCLUDE_TO_ADMINLOGINMIDDLEWARE', [])
        if view_name not in inclusion_set:
            return None

        if not request.user.is_superuser:
            return HttpResponseForbidden('You are not permitted')

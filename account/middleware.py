from django.shortcuts import render
from django.conf import settings

class IsEmailVarrifiedMiddleware(object):
    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request, *args, **kwargs):
        response = self.get_response(request)
        return response

    def process_view(self, request, view_func, *view_args, **view_kargs):
        view_name = '.'.join((view_func.__module__, view_func.__name__))
        # print(view_name)
        # If the view name is in our exclusion list, exit early
        inclusion_set = getattr(settings, 'INCLUDE_TO_IS_ACTIVE_MIDDLEWARE', [])
        if view_name not in inclusion_set:
            return None

        if not request.user.is_active:
            return render(request, 'account/email_verification.html')
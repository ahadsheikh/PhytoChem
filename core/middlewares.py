from django.http import HttpResponseNotFound


class AdminLoginMiddleware(object):
    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request, *args, **kwargs):
        print(request.user.is_superuser)
        if not request.user.is_superuser:
            return HttpResponseNotFound('You are not permitted')

        response = self.get_response(request)
        return response


def admin_login_middleware(get_response):

    def middleware(request):
        print(request.user.is_superuser)
        if not request.user.is_superuser:
            return HttpResponseNotFound('You are not permitted')

        response = get_response(request)
        return response

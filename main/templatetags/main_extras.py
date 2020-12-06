from django import template

register = template.Library()


def short_id(value):
    v = int(value.split('_')[1])
    return 'P' + str(v)


register.filter('short_id', short_id)


from django import template

register = template.Library()


# Make id short
def short_id(value):
    v = int(value.split('_')[1])
    return 'P' + str(v)


register.filter('short_id', short_id)


# Make a formula name as like real formula structure with subscript
def chem_formula(str):
    s_formula = ''
    for s in str:
        if s.isnumeric():
            s_formula = s_formula + '<sub>' + s + '</sub>'
        else:
            s_formula = s_formula + s
    return s_formula


register.filter('chem_formula', chem_formula)


def color_extract(arr, code):
    print(code)
    return arr[code];


register.filter('color_extract', color_extract)
from textwrap import indent


def indent_line(s, indent_with):
    if type(s) == str:
        return indent(s, indent_with, lambda line: True)
    else:
        return ''


def snake2camel(s):
    return ''.join(word.title() for word in s.split('_'))


def snake2pascal(s):
    camel = snake2camel(s)
    return camel[0].lower() + camel[1:]

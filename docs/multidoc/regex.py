
"""

Attributes
==========



p_cpp_tag : :func:`re.compile`
    .. code-block:: python

        p_cpp_tag = re.compile(
            r"(?P<indent>[^\S\\r\\n]+)?"   # Matches named group indent from start of line until comment.
            r"\/\/!"                     # Matches //! as user designed indentation of docstring.
            r"\s+"                       # Allows any number of whitespaces between start of comment and @.
            r"@get_docstring\("          # Matches '@get_docstring(' exactly.
            r"(\s+)?"                    # Allows any number (zero include) of whitespaces from '( until 'class' or 'func' name.
            r"("                         # Start of (optional) class regex group.
            r"(?P<cls>[\w]+)"                # Match a named group 'class' for a word of any length.
            r"\."                            # . match for separation between class name and method name.
            r"(?P<method>[\w]+)"             # Match a named group 'method' for a word of any length.
            r")?"                        # End of class regex group.
            r"(?P<func>\w+)?"                # Match (optional) named group 'func' for a word of any length.
            r"("                         # Start of overload matching group.
            r"(\s+)?"                        # Any number of whitespaces before comma.
            r","                             # Matches comma separating cls.method|func from overload variable.
            r"(\s+)?"                        # Any number (optional) of whitespaces after comma.
            r"(overload(\s+)?=(\s+)?"    # matches optional format of overload with overload key "overload = n"
            r")?"
            r"(?P<variant>[0-9]+)"       # Named group 'variant', int value between 0 and 9 ([0-9]), unlimited times (+)
            r")?"
            r"(\s+)?"                    # Optional whitespaces between closing parenthesis.
            r"\)",
            flags=re.MULTILINE)

    .. tip::

        See https://regex101.com/r/XJNQ8S/1 for desired effect.

p_python_tag : :func:`re.compile`

    .. code-block:: python

        p_python_tag = re.compile(
            r"(?P<indent>[^\S\\r\\n]+)?"   # Matches named group indent from start of line until comment.
            r"#"                         # Matches # as user designed indentation of docstring.
            r"\s+"                       # Allows any number of whitespaces between start of comment and @.
            r"@get_docstring\("          # Matches '@get_docstring(' exactly.
            r"(\s+)?"                    # Allows any number (zero include) of whitespaces from '( until 'class' or 'func' name.
            r"("                         # Start of (optional) class regex group.
            r"(?P<cls>[\w]+)"                # Match a named group 'class' for a word of any length.
            r"\."                            # . match for separation between class name and method name.
            r"(?P<method>[\w]+)"             # Match a named group 'method' for a word of any length.
            r")?"                        # End of class regex group.
            r"(?P<func>\w+)?"                # Match (optional) named group 'func' for a word of any length.
            r"("                         # Start of overload matching group.
            r"(\s+)?"                        # Any number of whitespaces before comma.
            r","                             # Matches comma separating cls.method|func from overload variable.
            r"(\s+)?"                        # Any number (optional) of whitespaces after comma.
            r"(overload(\s+)?=(\s+)?"    # matches optional format of overload with overload key "overload = n"
            r")?"
            r"(?P<variant>[0-9]+)"       # Named group 'variant', int value between 0 and 9 ([0-9]), unlimited times (+)
            r")?"
            r"(\s+)?"                    # Optional whitespaces between closing parenthesis.
            r"\)",
            flags=re.MULTILINE)

p_api_tag : :func:`re.compile`

    .. code-block:: python

        p_api_tag = re.compile(r".*#\s*\[(?P<expr>.*)\]")

p_package_file : :func:`re.compile`

    .. code-block:: python

        p_api_tag = re.compile(r".*__package__(.yml|.yaml)")

p_module_file : :func:`re.compile`

    .. code-block:: python

        p_module_file = re.compile(r".*(?P<module>\w+)(.yml|.yaml)")

"""
import re

p_cpp_tag = re.compile(
    r"(?P<indent>[^\S\r\n]+)?"   # Matches named group indent from start of line until comment.
    r"\/\/!"                     # Matches //! as user designed indentation of docstring.
    r"\s+"                       # Allows any number of whitespaces between start of comment and @.
    r"@get_docstring\("          # Matches '@get_docstring(' exactly.
    r"(\s+)?"                    # Allows any number (zero include) of whitespaces from '( until 'class' or 'func' name.
    r"("                         # Start of (optional) class regex group.
    r"(?P<cls>[\w]+)"                # Match a named group 'class' for a word of any length.
    r"\."                            # . match for separation between class name and method name.
    r"(?P<method>[\w]+)"             # Match a named group 'method' for a word of any length.
    r")?"                        # End of class regex group.
    r"(?P<func>\w+)?"                # Match (optional) named group 'func' for a word of any length.
    r"("                         # Start of overload matching group.
    r"(\s+)?"                        # Any number of whitespaces before comma.
    r","                             # Matches comma separating cls.method|func from overload variable.
    r"(\s+)?"                        # Any number (optional) of whitespaces after comma.
    r"(overload(\s+)?=(\s+)?"    # matches optional format of overload with overload key "overload = n"
    r")?"
    r"(?P<variant>[0-9]+)"       # Named group 'variant', int value between 0 and 9 ([0-9]), unlimited times (+)
    r")?"
    r"(\s+)?"                    # Optional whitespaces between closing parenthesis.
    r"\)",
    flags=re.MULTILINE)


p_python_tag = re.compile(
    r"(?P<indent>[^\S\r\n]+)?"   # Matches named group indent from start of line until comment.
    r"#"                         # Matches # as user designed indentation of docstring.
    r"\s+"                       # Allows any number of whitespaces between start of comment and @.
    r"@get_docstring\("          # Matches '@get_docstring(' exactly.
    r"(\s+)?"                    # Allows any number (zero include) of whitespaces from '( until 'class' or 'func' name.
    r"("                         # Start of (optional) class regex group.
    r"(?P<cls>[\w]+)"                # Match a named group 'class' for a word of any length.
    r"\."                            # . match for separation between class name and method name.
    r"(?P<method>[\w]+)"             # Match a named group 'method' for a word of any length.
    r")?"                        # End of class regex group.
    r"(?P<func>\w+)?"                # Match (optional) named group 'func' for a word of any length.
    r"("                         # Start of overload matching group.
    r"(\s+)?"                        # Any number of whitespaces before comma.
    r","                             # Matches comma separating cls.method|func from overload variable.
    r"(\s+)?"                        # Any number (optional) of whitespaces after comma.
    r"(overload(\s+)?=(\s+)?"    # matches optional format of overload with overload key "overload = n"
    r")?"                      
    r"(?P<variant>[0-9]+)"       # Named group 'variant', int value between 0 and 9 ([0-9]), unlimited times (+)
    r")?"
    r"(\s+)?"                    # Optional whitespaces between closing parenthesis.
    r"\)",
    flags=re.MULTILINE)


p_api_tag = re.compile(r".*#\s*\[(?P<expr>.*)\]")
p_package_file = re.compile(r".*__package__(.yml|.yaml)")
p_module_file = re.compile(r".*(?P<module>\w+)(.yml|.yaml)")


from document import *

import os

# clone repository
docstring_path = get_docstrings('https://github.com/tudat-team/tudat-multidoc.git')

# parse api declaration
api_declaration = parse_api_declaration(docstring_path, py=True)

# generate docstring header
docstrings_header_path = "/".join(os.path.abspath(__file__).split("/")[:-3])+"/tudatpy/kernel/docstrings.h"
generate_docstring_header(api_declaration, docstrings_header_path)
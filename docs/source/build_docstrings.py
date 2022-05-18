from document import *

import os

docstrings_path = "/".join(os.path.abspath(__file__).split("/")[:-3])+"/include/tudatpy/docstrings.h"

# clone repository
docstring_path = get_docstrings('https://github.com/tudat-team/tudat-multidoc.git')

# parse api declaration
api_declaration = parse_api_declaration(docstring_path, py=True)

# generate docstring header
generate_docstring_header(api_declaration, docstrings_path)
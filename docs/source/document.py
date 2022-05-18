# insert current directory into path
import sys, os
sys.path.insert(0, "/".join(os.path.abspath(__file__).split("/")[:-2]))

import tempfile
from multidoc.generate import generate_py_sphinx, generate_docstring_header
from multidoc.parsing import parse_api_declaration


# get docstrings from github
def get_docstrings(git_url, git_rev=None):
    import os
    import subprocess
    import shutil
    import tempfile

    # create temporary directory
    tmp_dir = tempfile.mkdtemp()

    # clone repository
    subprocess.check_call(['git', 'clone', git_url, tmp_dir])

    # checkout revision
    if git_rev is not None:
        subprocess.check_call(['git', 'checkout', git_rev], cwd=tmp_dir)

    # delete docstrings folder if exists
    if os.path.exists('docstrings'):
        shutil.rmtree('docstrings')

    # copy files to destination
    shutil.copytree(os.path.join(tmp_dir, 'docstrings'), 'docstrings')
    # shutil.copytree(os.path.join(tmp_dir, 'multidoc'), '../multidoc')

    # remove temporary directory
    shutil.rmtree(tmp_dir)

    # return absolute docstring path
    return os.path.abspath('docstrings')

# generate documentation in temp directory
def generate_documentation(api_declaration, output_dir):
    import os
    import shutil

    # create temporary directory
    tmp_dir = tempfile.mkdtemp()

    # copy files to destination
    generate_py_sphinx(api_declaration=api_declaration,
                       dest_dir=tmp_dir)

    # copy files to destination
    shutil.copytree(os.path.join(tmp_dir), output_dir, dirs_exist_ok=True)

    # remove temporary directory
    shutil.rmtree(tmp_dir)

    # return documentation source path
    return os.path.abspath(output_dir)


if __name__ == '__main__':
    multidoc_git_url = 'https://github.com/tudat-team/tudat-multidoc.git'
    multidoc_git_rev = '0811926d9f98331a5d0eca0108e44e7acb6a972c'

    # clone repository
    docstring_path = get_docstrings(multidoc_git_url, multidoc_git_rev)

    # parse api declaration
    api_declaration = parse_api_declaration(docstring_path, py=True)

    # source path
    source_path = generate_documentation(api_declaration, '.')
import shutil
import os
import re
import jinja2
from pathlib import Path
from multidoc.template import TEMPLATE_DIR, TEMPLATE_DOCSTRING_HEADER
from multidoc.parsing import parse_api_declaration
from multidoc.regex import p_cpp_tag
from multidoc.utils import parts, indent_line, snake2pascal, recurse_dict
from multidoc.error import *


def _generate_documented(target_src, force_overwrite=True, dest=None):
    path = Path(target_src)
    base = os.path.basename(path)
    if dest is None:
        dest = os.path.join(path, "..", "." + base + "-documented")
    if force_overwrite:
        if os.path.exists(dest):
            shutil.rmtree(dest)
    shutil.copytree(target_src, dest)
    return dest, base


def generate_docstring_header(api_declaration,
                              destination,
                              template_path=TEMPLATE_DOCSTRING_HEADER):
    """

    Parameters
    ----------
    api_declaration : APIDeclaration
    destination
    template_path

    Returns
    -------

    """
    # read docstrings template
    with open(template_path, "r") as f:
        t = jinja2.Template(f.read(),
                            lstrip_blocks=True,
                            trim_blocks=True)

    # generate docstring file
    with open(destination, "w") as f:
        import json
        with open("json.json", "w") as ff:
            ff.write(json.dumps(api_declaration, indent=4))
        f.write(t.render(api_structure=api_declaration))


def generate_cpp_docstring(api_declaration, include_path, dest):
    # load structure from api definition
    structure = api_declaration

    for path in Path(include_path).rglob("*.h"):
        _parts = parts(path)
        include_idx = _parts.index(os.path.basename(include_path))
        module_structure = _parts[include_idx + 1:]
        assert module_structure[0] == structure[
            "name"]  # currently, include/<package_name> is mandatory

        # i.e. a, *b, c = [0, 1, 2, 3, 4] -> a=0; b=[1,2,3]; c=4
        package_name, *modules, filename = module_structure

        # this is a very temporary solution to dealing with looking for
        #  docstrings that do not exist.
        try:
            t_structure = recurse_dict(structure, modules)
        except KeyError:
            continue

        # open the .hpp file and replace all tags with docstrings - if found.
        with open(path, "r") as f:
            processed = raw = f.read()
            for match in re.finditer(p_cpp_tag, raw):
                c = match.group('cls') if match.group('cls') else None
                m = match.group('method') if match.group('method') else None
                v = match.group('variant') if match.group('variant') else None
                f = match.group('func') if match.group('func') else None
                i = match.group('indent') if match.group('indent') else ''
                if c:
                    try:
                        aux = t_structure[c]
                    except KeyError:
                        raise ClassNotDeclaredError(f"{c} was not found in "
                                                    f"{'/'.join(modules)} for "
                                                    f"{package_name} API "
                                                    f"declaration.")
                    try:
                        aux = aux[m]
                    except KeyError:
                        raise MethodNotDeclaredError(f"{m} was not found in "
                                                     f"{'/'.join(modules + [c])} for "
                                                     f"{package_name} API "
                                                     f"declaration.")
                elif f:
                    try:
                        aux = t_structure[f]
                    except KeyError:
                        raise FunctionNotDeclaredError(f"{f} was not found in "
                                                       f"{'/'.join(modules)} for "
                                                       f"{package_name} API ")
                if v:
                    try:
                        aux = aux[int(v)]
                    except KeyError:
                        raise OverloadNotFoundError(
                            f"Overload ({v}) for {aux['name']} was "
                            f"not found, are the overloads stored as"
                            f" lists, as they should be?")
                processed = processed.replace(match.group(0),
                                              indent_line(aux, i))

        # create destination directory structure, if it does not exist.
        if not os.path.exists(os.path.join(dest, *modules)):
            os.makedirs(os.path.join(dest, *modules))

        # save the .hpp file with all tags replaced with docstrings.
        with open(os.path.join(dest, *modules, filename), "w") as f2:
            f2.write(processed)


def generate_cpp_documented(api_prefix, target_src, dest=None):
    path_documented, name = _generate_documented(target_src, dest=dest)
    api_declaration=parse_api_declaration(api_prefix, cpp=True)
    generate_cpp_docstring(api_declaration=api_declaration,
                           dest=os.path.join(path_documented,
                                             f"include/{name}"),
                           include_path=os.path.join(target_src, "include"))
    generate_cpp_sphinx(api_declaration=api_declaration,
                        dest_dir=os.path.join(path_documented, 'docs',
                                              'sphinx', 'source'))


def generate_cpp_sphinx(
        dest_dir, api_declaration,
        template=os.path.join(TEMPLATE_DIR, "cpp-sphinx-module.jinja2")):
    if not os.path.exists(dest_dir):
        # shutil.rmtree(dest_dir)
        os.makedirs(dest_dir)

    # load api declaration structure
    structure = api_declaration
    # print(json.dumps(structure, indent=4))  # for inspection of structure
    with open(template, "r") as f:
        t = jinja2.Template(f.read())

    t.globals["repeat_str"] = lambda n, s: len(n) * s

    namespace_list = [structure['name']]

    # Generate index.rst
    s = t.render(structure=structure, namespace_list=namespace_list,
                 title="API Reference")
    t.globals["snake2pascal"] = snake2pascal
    with open(os.path.join(dest_dir, f"index.rst"), "w") as f:
        f.write(s)

    def recurse(structure, namespace_list):
        title = '``' + "::".join(namespace_list) + '``'
        if "modules" in structure.keys():
            for module in structure['modules']:
                recurse(structure[module], namespace_list + [module])
            s = t.render(structure=structure, namespace_list=namespace_list,
                         title=title)
            with open(os.path.join(dest_dir, f"{structure['name']}.rst"),
                      "w") as f:
                f.write(s)

        else:
            s = t.render(structure=structure, namespace_list=namespace_list,
                         title=title)
            with open(os.path.join(dest_dir, f"{structure['name']}.rst"),
                      "w") as f:
                f.write(s)

    for module in structure["modules"]:
        recurse(structure[module], namespace_list=namespace_list + [module])
    #


def generate_py_sphinx(
        dest_dir, api_declaration,
        template=os.path.join(TEMPLATE_DIR, "py-sphinx-module.jinja2")):
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)

    # load api declaration structure
    structure = api_declaration

    # read pybind docstrings template
    with open(template, "r") as f:
        t = jinja2.Template(f.read())

    t.globals["repeat_str"] = lambda n, s: len(n) * s

    namespace_list = [structure['name']]

    # Generate index.rst
    s = t.render(structure=structure, namespace_list=namespace_list,
                 title="API Reference")
    with open(os.path.join(dest_dir, f"index.rst"), "w") as f:
        f.write(s)

    def recurse(structure, namespace_list):
        title = '``' + ".".join(namespace_list[-1:])+ '``'
        if "modules" in structure.keys():
            for module in structure['modules']:
                recurse(structure[module], namespace_list + [module])
            s = t.render(structure=structure, namespace_list=namespace_list,title=title)
            with open(os.path.join(dest_dir, f"{structure['name']}.rst"),"w") as f:
                f.write(s)

        else:
            s = t.render(structure=structure, namespace_list=namespace_list,title=title)
            with open(os.path.join(dest_dir, f"{structure['name']}.rst"),"w") as f:
                f.write(s)

    for module in structure["modules"]:
        recurse(structure[module], namespace_list=namespace_list + [module])


def generate_pybind_documented(api_prefix, target_src, dest=None):
    api_declaration = parse_api_declaration(api_prefix, py=True)
    path_documented, name = _generate_documented(target_src, dest=dest)
    generated_header_dir = os.path.join(path_documented, f"include/{name}")
    if not os.path.exists(generated_header_dir):
        os.makedirs(generated_header_dir)
    generate_docstring_header(api_declaration, os.path.join(generated_header_dir,
                                                       "docstrings.h"))
    # print(os.path.join(path_documented, 'docs', 'source'))
    generate_py_sphinx(api_declaration=api_declaration,
                       dest_dir=os.path.join(path_documented, 'docs',
                                             'source'))


if __name__ == "__main__":
    generate_cpp_sphinx(
        api_prefix="../tests/test-docstrings",
        dest_dir="test-sphinx-cpp")

    generate_py_sphinx(
        api_prefix="../tests/test-docstrings",
        dest_dir="test-sphinx-py")

import re
import os
from multidoc.template import render_python_docstring
from multidoc.template import render_cpp_docstring
from multidoc.template import get_docstring_template
from multidoc.parsing.io import yaml2dict
from multidoc.parsing import logger
from multidoc.regex import p_package_file, p_module_file
import pathlib

def guess_project_type(project_src):
    logger.warning(f"GUESSING project type, because `project_type` was not found for src {project_src}!")
    files_cpp = pathlib.Path(project_src).glob("*.cpp|*.h|*.hpp")
    files_py = pathlib.Path(project_src).glob("*.py")

    # should check cpp for pybind headers.
    # header_pybind =
    if files_cpp and files_py:
        return "py"
    elif files_cpp:
        return "cpp"
    elif files_py:
        return "py"


def parse_function(function, local):
    """

    Parameters
    ----------
    function
    local

    Returns
    -------

    """
    logger.info(f"Parsing function: {function.name} with locals: {local}")
    if "cpp" in local.keys() and "py" in local.keys():
        assert (local["cpp"] == local["cpp"] & local[
            "cpp"] is True) is False
    if "cpp" in local.keys():
        if local["cpp"]:
            return render_cpp_docstring(function)

    if "py" in local.keys():
        if local["py"]:
            return render_python_docstring(function)


from collections import defaultdict


def parse_functions(structure, functions, **kwargs):
    nl = '\n'
    logger.info(
        f"Parsing the following functions with {kwargs}: {nl}{nl.join([r'    - '+ function.name for function in functions] + [''])} "
    )
    t = get_docstring_template(**kwargs)
    result = defaultdict(list)
    for i, f in enumerate(functions):
        result[f.name].append(i)

    for k, v in result.items():
        if len(v) > 1:  # overloaded
            structure[k] = defaultdict(list)
            structure[k]["overloaded"] = True
            structure[k]["overloads"] = len(v)
            for idx, overload in enumerate(v):
                structure[k][idx] = t.render(**functions[overload].dict())
            structure[k] = dict(structure[k])
        else:
            structure[k] = t.render(**functions[v[0]].dict())

    return structure

    # logger.info(f"Parsing class: {cls.name} with locals: {local}")
    # _return = {}
    # for method in cls.methods:
    #     _return[method.name] = parse_method(method, local)
    # _return["__docstring__"] = parse_method(cls, local)
    # return _return


def parse_classes(structure, classes, **kwargs):
    nl = '\n'
    logger.info(
        f"Parsing the following classes with {kwargs}: {nl}{nl.join([r'    - '+cls.name for cls in classes] + [''])} "
    )
    t = get_docstring_template(**kwargs)
    for cls in classes:
        _result = {}
        if cls.methods:
            _result.update(parse_functions(cls.dict(), cls.methods, **kwargs))
        _result.update({"__docstring__": t.render(**cls.dict())})
        structure[cls.name] = _result

    return structure


def parse_method(function, local):
    """

    Parameters
    ----------
    function
    local

    Returns
    -------

    """
    logger.info(f"Parsing method: {function.name} with locals: {local}")
    if "cpp" in local.keys() and "py" in local.keys():
        assert (local["cpp"] == local["cpp"] & local[
            "cpp"] is True) is False
    if "cpp" in local.keys():
        if local["cpp"]:
            return render_cpp_docstring(function)

    if "py" in local.keys():
        if local["py"]:
            return render_python_docstring(function)


import json


def parse_class(cls, local):
    """

    Parameters
    ----------
    cls
    local

    Returns
    -------

    """
    logger.info(f"Parsing class: {cls.name} with locals: {local}")
    _return = {}
    for method in cls.methods:
        _return[method.name] = parse_method(method, local)
    _return["__docstring__"] = parse_method(cls, local)
    return _return


def parse_constant(constant, local):
    """

    Parameters
    ----------
    constant
    local

    Returns
    -------

    """
    return constant


from multidoc.parsing.models import Module, Package


def parse_api_declaration(path: str, parent=None, **kwargs):
    """

    Parameters
    ----------
    path
    local
    parent

    Returns
    -------

    """
    structure = dict()
    local = kwargs if kwargs is not None else dict()

    # if directory given, treats it as package declaration
    if os.path.isdir(path):

        # look for package declaration __package__(.yml|.yaml)
        files = list(filter(p_package_file.match, os.listdir(path)))

        if len(files) == 0:
            raise ModuleNotFoundError("__package__.yaml/yml not found in "
                                      "directory path.")
        elif len(files) > 1:
            raise ModuleNotFoundError("Multiple __package__.yaml/yml files "
                                      "found in directory path.")
        else:
            module = Package.parse_yaml(os.path.join(path, files[0]), **kwargs)

            # define type as package.
            structure["type"] = "package"
            structure["path"] = path
            structure["file"] = files[0]
            structure["_implicit_name"] = os.path.basename(path)

    # if file given, treats it as module declaration
    elif os.path.isfile(path):

        _, extension = os.path.splitext(path)
        #
        if extension == ".yaml" or extension == ".yml":

            if not p_package_file.match(path):
                module = Module.parse_yaml(path)

                # define type as package.
                structure["type"] = "module"
                structure["path"] = os.path.dirname(path)
                structure["file"] = os.path.basename(path)
                structure["_implicit_name"] = \
                    os.path.split(os.path.basename(path))[0]
            else:
                module = Package.parse_yaml(path)
                structure["type"] = "package"
                structure["path"] = os.path.dirname(path)
                structure["file"] = os.path.basename(path)
                structure["_implicit_name"] = \
                    os.path.split(os.path.dirname(path))[0]
        else:
            raise ModuleNotFoundError("Only .yml or .yaml files can be used "
                                      "to declare modules and packages.")

    # then path was given without yml/yaml extension.
    # (i.e. directory/module, where directory/module.yml/yaml exists)
    # TODO: Is this really necessary?
    elif list(
            filter(re.compile(fr"{os.path.basename(path)}(.yml|.yaml)").match,
                   os.listdir(os.path.dirname(path)))):

        matches = list(
            filter(re.compile(fr"{os.path.basename(path)}(.yml|.yaml)").match,
                   os.listdir(os.path.dirname(path))))
        basename = os.path.basename(path)

        if len(matches) > 1:
            raise ModuleNotFoundError(f"Multiple {basename}.yaml/yml files "
                                      "found for given path.")
        else:
            module = Module.parse_yaml(
                os.path.join(os.path.dirname(path), matches[0]), **kwargs)

            # define type as package.
            structure["type"] = "module"
            structure["path"] = os.path.dirname(path)
            structure["file"] = matches[0]
            structure["_implicit_name"] = os.path.split(
                os.path.basename(path))[-1]

    # user gives path that is neither file nor directory.
    else:
        raise ModuleNotFoundError("Path provided was not recognized as a "
                                  "directory containing a __package__.yaml/yml"
                                  "or as a module.yml/.yaml declaration.")

    structure.update(module.dict())

    if module.config:
        # parse multidoc related configuration
        if module.config.name:  # package
            structure["name"] = module.config.name
        else:
            structure["name"] = structure["_implicit_name"]

        if module.config.version:
            structure["version"] = module.config.version
        else:
            structure["version"] = None
    else:
        structure["name"] = structure["_implicit_name"]
        structure["version"] = None

    # module level functions
    if module.functions:
        structure = parse_functions(structure, module.functions, **local)

    # module level functions
    if module.classes:
        structure = parse_classes(structure, module.classes, **local)

    # module level functions
    if module.constants:
        for constant in module.constants:
            structure[constant.name] = parse_constant(constant, local)

    # iterate through api level modules
    if type(module) == Package:
        for submodule in module.modules:
            module_path = os.path.join(structure["path"], submodule)
            structure[submodule] = parse_api_declaration(module_path,
                                                         parent=structure,
                                                         **local)
            # TODO: Figure out how to deal with submodule name clashes with
            #  reserved names. Hide all docstrings behind __docs__ key?

    return structure


if __name__ == "__main__":
    s = parse_api_declaration("../../tests/test-docstrings",
                              py=True)
    if s["summary"]:
        print(s["summary"])
    # print(s)
    import json

    json.dumps(s, indent=4)

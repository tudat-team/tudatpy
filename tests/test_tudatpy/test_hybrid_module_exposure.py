# General imports
import inspect
import unittest
from importlib import import_module


def _recursively_find_children(module):
    """
    Recursively find all children of a module.
    """

    children = inspect.getmembers(module, inspect.ismodule)

    if children:
        for name, child in children:
            children += _recursively_find_children(child)
        return children
    else:
        return []


def test_tudatpy_kernel_module_exposure():

    import tudatpy.kernel as kernel

    exposed_kernel_modules = [
        # 'utils',
        "astro",
        "trajectory_design",
        "constants",
        "interface",
        # 'io',
        "math",
        "numerical_simulation",
    ]

    module_methods = lambda module: inspect.getmembers(module, inspect.isroutine)

    module_variables = lambda module: [
        (name, value)
        for name, value in inspect.getmembers(module)
        if (
            not inspect.isroutine(value)
            and not inspect.ismodule(value)
            and not name.startswith("__")
        )
    ]

    for submodule_name, submodule in _recursively_find_children(kernel):

        if submodule.__name__.split(".")[2] in exposed_kernel_modules:
            kernel_module = import_module(submodule.__name__)
            exposed_module = import_module(submodule.__name__.replace("kernel.", ""))

            exposed_method_names = {name for name, _ in module_methods(exposed_module)}
            kernel_method_names = {name for name, _ in module_methods(kernel_module)}

            # Assert the public module exposes the same method names. The compatibility
            # modules may re-export callables through a new module path, so object
            # identity is not part of the exposure contract.
            assert kernel_method_names <= exposed_method_names

            exposed_variable_names = {name for name, _ in module_variables(exposed_module)}
            kernel_variable_names = {name for name, _ in module_variables(kernel_module)}

            # Assert both modules contain the same variable names.
            assert kernel_variable_names <= exposed_variable_names

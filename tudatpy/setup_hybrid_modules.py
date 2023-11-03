''' 
Copyright (c) 2010-2023, Delft University of Technology
All rigths reserved

This file is part of the Tudat. Redistribution and use in source and 
binary forms, with or without modification, are permitted exclusively
under the terms of the Modified BSD license. You should have received
a copy of the license with this file. If not, please or visit:
http://tudat.tudelft.nl/LICENSE.
'''

# General imports
import os
import re
import inspect
import argparse
from pathlib import Path

#--------------------------------------------------------------------
#%% FUNCTION TO EXPOSE ALL HYBRID KERNEL MODULES
#--------------------------------------------------------------------

disclaimer = lambda module_name: f"""# THE FOLLOWING IMPORT STATEMENT HAS BEEN AUTOMATICALLY GENERATED, DO NOT EDIT DIRECTLY
#
# This file, by virtue of the import statement below, merges
# the Tudat kernel module `tudatpy.kernel.{module_name}` with
# its Python extensions defined in `tudatpy/{module_name.replace('.', '/')}`.
# 
# This allows the import of all the C++ and Python submodules of the 
# `{module_name}` kernel module directly from tudatpy:
# 
#     from tudatpy.{module_name} import <any>
# 
# Without the statement below, importing the `{module_name}` kernel module
# would only be possible as follows, and hybrid Python/C++ modules would not
# be posible in tudatpy.
# 
#     from tudatpy.kernel.{module_name} import <any>
# 
# The reason why C++ kernel modules can only be imported as written above 
# is an issue with the `def_submodule` function of pybind11. The issue is discussed
# [here](https://github.com/pybind/pybind11/issues/2639). 
# 
# We circumvent the issue by automatically creating an *empty* Python module for each
# kernel module and submodule, and exposing the kernel module or submodule from the
# Python module by adding the import statement below to the Python module's `__init__.py` (this file). 
# This workaround was proposed in [this comment](https://github.com/pybind/pybind11/issues/2639#issuecomment-721238757).
# 
# An added benefit of this method is that it makes it possible to write Python extensions 
# and add them to the kernel modules simply by placing them inside the newly created Python 
# module (at the same level as this file), which is not possible with the `def_submodule` 
# function of pybind11.
# This is automatically done at compile-time, by moving the Python extensions of this module 
# from `tudatpy/{module_name.replace('.', '/')}` to this module (the module
# corresponding to this `__init__.py` file: `tudatpy.{module_name}`). 
"""

def expose_hybrid_module(module):
    """
    Expose all Tudat kernel modules as hybrid C++/Python tudatpy modules.
    """

    module_path  = module.__name__
    module_name  = module_path[len('tudatpy.kernel.'):]
    module_depth = max(len(module_name.split('.')) - 3, 0)

    # Hybrid module path
    hybrid_module_path = Path(f"tudatpy/{module_name.replace('.', '/')}")
    hybrid_module_path.mkdir(parents=True, exist_ok=True)

    # __init__.py
    import_statement = f"from {module_path} import *"
    disclaimer_text  = disclaimer(module_name)
    init_file_path   = hybrid_module_path / '__init__.py'

    if not os.path.isfile(init_file_path):
        # If __init__.py does not exist, create it directly
        with open(init_file_path, 'w') as f:
            f.write(f"{disclaimer_text}\n{import_statement}")
    else:
        # Otherwise, if the import statement is not already present in __init__.py,
        # append the disclaimer and import statement to __init__.py
        with open(init_file_path, 'r') as f:
            contents = f.read()
        if re.search(re.escape(import_statement), contents):
            # import_statement is already present in __init__.py
            pass
        else:
            with open(init_file_path, 'a') as f:
                f.write(f"\n\n{disclaimer_text}\n{import_statement}")

#--------------------------------------------------------------------
#%% OBTAIN LIST CONTAINING ALL KERNEL MODULES AND SUBMODULES
#--------------------------------------------------------------------

def recursively_find_children(module):
    """
    Recursively find all children of a module.
    """

    children = inspect.getmembers(module, inspect.ismodule)

    if children:
        for (name, child) in children:
            children += recursively_find_children(child)
        return children
    else:
        return []

#--------------------------------------------------------------------
#%% EXPOSE KERNEL MODULES AND SUBMODULES INPUT BY USER
#--------------------------------------------------------------------

if __name__ == '__main__':

    # Import Tudat kernel
    import tudatpy.kernel as kernel

    parser = argparse.ArgumentParser(description='Expose hybrid kernel modules.')
    parser.add_argument('modules', metavar='module', type=str, nargs='+',
                        help='a list of hybrid kernel modules to expose')
    parser.add_argument('-v', '--verbose', action=argparse.BooleanOptionalAction,
                        help='whether to announce the module and submodule being exposed')
    args = parser.parse_args()

    for module_name in args.modules:

        if args.verbose:
            message = f"EXPOSING MODULE: {module_name}"
            message_width = 100
            print(f"\n{' ' + message + ' ':=^100}")

        # Identify and expose module
        module = getattr(kernel, module_name)
        expose_hybrid_module(module)
        # Find and expose all submodules
        for (submodule_name, submodule) in recursively_find_children(module):
            # Expose module
            expose_hybrid_module(submodule)
            # Report completion
            if args.verbose:
                print(f"{'-':>40} {submodule_name}")
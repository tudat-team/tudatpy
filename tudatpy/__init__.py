'''
Copyright (c) 2010-2024, Delft University of Technology
All rigths reserved

This file is part of the Tudat. Redistribution and use in source and
binary forms, with or without modification, are permitted exclusively
under the terms of the Modified BSD license. You should have received
a copy of the license with this file. If not, please or visit:
http://tudat.tudelft.nl/LICENSE.
'''

# Stdlib imports
import os
import sys
import glob
import importlib
import importlib.abc
import importlib.util
import importlib.machinery
from setuptools import find_packages

# Import tudat kernel
from tudatpy import kernel

# Tudatpy version information
from tudatpy._version import *

# Import Tudat version information
from tudatpy.kernel import \
    _tudat_version, \
    _tudat_version_major, \
    _tudat_version_minor, \
    _tudat_version_patch

"""
Expose all tudat kernel modules directly from tudatpy
"""

# List all tudatpy Python packages
tudatpy_path = os.path.dirname(__file__)
tudatpy_packages = find_packages(tudatpy_path)
tudatpy_packages = [f'tudatpy.{module}' for module in tudatpy_packages]
tudatpy_modules = ['tudatpy.' + f[len(tudatpy_path)+1:-len('.py')].replace('/', '.') for f in glob.glob(os.path.join(tudatpy_path, '**/*.py'), recursive=True) if '__init__' not in f]

class TudatpyModuleFinder(importlib.abc.MetaPathFinder):
    def load_tudatpy_module(self, fullname, path, target):
        for finder in sys.meta_path:
            if finder is not self:
                tudatpy_module_spec = finder.find_spec(fullname, path, target)
                if tudatpy_module_spec is not None:
                    break
        tudatpy_module = importlib.util.module_from_spec(tudatpy_module_spec)
        tudatpy_module_spec.loader.exec_module(tudatpy_module)
        return tudatpy_module
    def find_spec(self, fullname, path, target=None):
        if fullname.startswith("tudatpy"):

            # Determine whether module is a Tudatpy Python module
            is_tudatpy_module = fullname in tudatpy_packages + tudatpy_modules

            # Determine whether module is a Tudat kernel module
            kernel_module_name = fullname.replace('tudatpy', 'tudatpy.kernel', 1)
            try:
                kernel_module = importlib.import_module(kernel_module_name)
            except Exception as e:
                kernel_module = None
            is_kernel_module = kernel_module is not None

            # If the module is a kernel module, wrap it
            if is_kernel_module and not is_tudatpy_module:
                spec = importlib.util.spec_from_loader(
                    name=fullname,
                    loader=KernelModuleWrapper(kernel_module)
                )
                return spec
            # If the module is a hybrid module
            elif is_kernel_module and is_tudatpy_module:
                # Create default spec
                for finder in sys.meta_path:
                    if finder is not self:
                        tudatpy_module_spec = finder.find_spec(fullname, path, target)
                        if tudatpy_module_spec is not None:
                            break
                # Return hybrid module spec
                hybrid_module_spec = importlib.util.spec_from_loader(
                    name=fullname,
                    loader=HybridModuleLoader(
                        tudatpy_module_spec,
                        kernel_module
                    )
                )
                return hybrid_module_spec
            else:
                return None
        return None

class KernelModuleWrapper(importlib.abc.Loader):
    def __init__(self, module) -> None:
        self.kernel_module = module
    def create_module(self, spec: importlib.machinery.ModuleSpec) -> None:
        return self.kernel_module
    def exec_module(self, module) -> None:
        module.__path__ = __path__
        sys.modules[module.__name__] = self.kernel_module

class HybridModuleLoader(importlib.abc.Loader):
    def __init__(self, tudatpy_module_spec, kernel_module) -> None:
        self.tudatpy_module_spec = tudatpy_module_spec
        self.kernel_module = kernel_module
    def create_module(self, spec: importlib.machinery.ModuleSpec) -> None:
        tudatpy_module = importlib.util.module_from_spec(self.tudatpy_module_spec)
        # Add module to sys.modules
        sys.modules[tudatpy_module.__name__] = tudatpy_module
        # Execute tudatpy module
        self.tudatpy_module_spec.loader.exec_module(tudatpy_module)
        return tudatpy_module
    def exec_module(self, module) -> None:
        # Set module __path__
        module.__path__ = getattr(module, '__path__', [])
        # Expose kernel module
        def getattr_from_hybrid_module(attr: str):
            try:
                return self.kernel_module.__dict__[attr]
            except KeyError:
                raise AttributeError(f"neither module '{module.__name__}' nor Tudat kernel module '{module.__name__.replace('tudatpy', 'tudatpy.kernel')}' have no attribute {attr}")
        module.__getattr__ = getattr_from_hybrid_module

# Install the custom finder
sys.meta_path.insert(0, TudatpyModuleFinder())

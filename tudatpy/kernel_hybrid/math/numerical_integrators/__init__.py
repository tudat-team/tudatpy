# This file, by virtue of the import statement below, merges
# the Tudat kernel module `tudatpy.kernel.math.numerical_integrators` with
# its Python extensions defined in `tudatpy/kernel_hybrid/math/numerical_integrators`.
# 
# This allows the import of all the C++ and Python submodules of the 
# `math.numerical_integrators` kernel module directly from tudatpy:
# 
#     from tudatpy.math.numerical_integrators import <any>
# 
# Without the statement below, importing the `math.numerical_integrators` kernel module
# would only be possible as follows, and hybrid Python/C++ modules would not
# be posible in tudatpy.
# 
#     from tudatpy.kernel.math.numerical_integrators import <any>
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
# An added benefit of this method is that it makes it possible to write Python extensions 
# and add them to the kernel modules simply by placing them inside this module!

from tudatpy.kernel.math.numerical_integrators import *
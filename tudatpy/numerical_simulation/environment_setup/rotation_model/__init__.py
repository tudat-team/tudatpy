# This file, by virtue of the import statement below, merges
# the Tudat kernel module `tudatpy.kernel.numerical_simulation.environment_setup.rotation_model` with
# its Python extensions defined in `tudatpy/numerical_simulation/environment_setup/rotation_model`.
# 
# This allows the import of all the C++ and Python submodules of the 
# `numerical_simulation.environment_setup.rotation_model` kernel module directly from tudatpy:
# 
#     from tudatpy.numerical_simulation.environment_setup.rotation_model import <any>
# 
# Without the statement below, importing the `numerical_simulation.environment_setup.rotation_model` kernel module
# would only be possible as follows, and hybrid Python/C++ modules would not
# be posible in tudatpy.
# 
#     from tudatpy.kernel.numerical_simulation.environment_setup.rotation_model import <any>
# 
# The reason why C++ kernel modules can only be imported as written above 
# is an issue with the `def_submodule` function of pybind11. The issue is discussed
# [here](https://github.com/pybind/pybind11/issues/2639). 
# 
# We circumvent the issue by, for each module and submodule, 
# 
#     1. Creating a Python module (an empty directory with an `__init__.py` file) 
#        inside `tudatpy` with the same name
#     2. Adding the import statement below to the Python module's `__init__.py` 
#        (this file), thereby making the kernel module or submodule from the Python module.
# 
# This workaround was proposed in [this comment](https://github.com/pybind/pybind11/issues/2639#issuecomment-721238757).
# 
# An added benefit of this method is that it makes it possible to write Python extensions 
# and add them to the kernel modules simply by placing them inside the newly created Python 
# module (at the same level as this file), which is not possible with the `def_submodule` 
# function of pybind11.
# An added benefit of this method is that it makes it possible to write Python extensions 
# and add them to the kernel modules simply by placing them inside this module!

from tudatpy.kernel.numerical_simulation.environment_setup.rotation_model import *
from tudatpy.numerical_simulation.environment_setup.rotation_model._import_all_kernel_members import *
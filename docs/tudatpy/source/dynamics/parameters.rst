.. _parameters:

``parameters``
==============

Functions
---------
.. currentmodule:: tudatpy.dynamics.parameters

.. autosummary::

   print_parameter_names
   set_parameter_correction_constraint
   set_initial_state_correction_constraint
   remove_parameter_correction_constraint
   remove_parameter_correction_constraints

.. autofunction:: tudatpy.dynamics.parameters.print_parameter_names

.. autofunction:: tudatpy.dynamics.parameters.set_parameter_correction_constraint

.. autofunction:: tudatpy.dynamics.parameters.set_initial_state_correction_constraint

.. autofunction:: tudatpy.dynamics.parameters.remove_parameter_correction_constraint

.. autofunction:: tudatpy.dynamics.parameters.remove_parameter_correction_constraints



Classes
-------
.. currentmodule:: tudatpy.dynamics.parameters

.. autosummary::

   EstimatableParameter

   EstimatableParameterSet

.. autoclass:: tudatpy.dynamics.parameters.EstimatableParameter
   :members:

.. autosummary::

   EstimatableParameter.get_parameter_value
   EstimatableParameter.set_parameter_value

.. autoclass:: tudatpy.dynamics.parameters.EstimatableParameterSet
   :members:

.. autosummary::

   EstimatableParameterSet.indices_for_parameter_type
   EstimatableParameterSet.parameters_for_parameter_type

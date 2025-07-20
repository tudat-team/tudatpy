.. _integrator:

``integrator``
==============
This module provides the functionality for creating integrator
settings.












Functions
---------
.. currentmodule:: tudatpy.dynamics.propagation_setup.integrator

.. autosummary::

   runge_kutta_fixed_step

   runge_kutta_variable_step

   bulirsch_stoer_fixed_step

   bulirsch_stoer_variable_step

   adams_bashforth_moulton

   adams_bashforth_moulton_fixed_order

   adams_bashforth_moulton_fixed_step

   adams_bashforth_moulton_fixed_step_fixed_order

   print_butcher_tableau

   step_size_validation

   step_size_control_elementwise_scalar_tolerance

   step_size_control_elementwise_matrix_tolerance

   step_size_control_blockwise_scalar_tolerance

   step_size_control_blockwise_matrix_tolerance

   step_size_control_custom_blockwise_scalar_tolerance

   step_size_control_custom_blockwise_matrix_tolerance

   standard_cartesian_state_element_blocks





.. autofunction:: tudatpy.dynamics.propagation_setup.integrator.runge_kutta_fixed_step

.. autofunction:: tudatpy.dynamics.propagation_setup.integrator.runge_kutta_variable_step

.. autofunction:: tudatpy.dynamics.propagation_setup.integrator.bulirsch_stoer_fixed_step

.. autofunction:: tudatpy.dynamics.propagation_setup.integrator.bulirsch_stoer_variable_step

.. autofunction:: tudatpy.dynamics.propagation_setup.integrator.adams_bashforth_moulton

.. autofunction:: tudatpy.dynamics.propagation_setup.integrator.adams_bashforth_moulton_fixed_order

.. autofunction:: tudatpy.dynamics.propagation_setup.integrator.adams_bashforth_moulton_fixed_step

.. autofunction:: tudatpy.dynamics.propagation_setup.integrator.adams_bashforth_moulton_fixed_step_fixed_order

.. autofunction:: tudatpy.dynamics.propagation_setup.integrator.print_butcher_tableau

.. autofunction:: tudatpy.dynamics.propagation_setup.integrator.step_size_validation

.. autofunction:: tudatpy.dynamics.propagation_setup.integrator.step_size_control_elementwise_scalar_tolerance

.. autofunction:: tudatpy.dynamics.propagation_setup.integrator.step_size_control_elementwise_matrix_tolerance

.. autofunction:: tudatpy.dynamics.propagation_setup.integrator.step_size_control_blockwise_scalar_tolerance

.. autofunction:: tudatpy.dynamics.propagation_setup.integrator.step_size_control_blockwise_matrix_tolerance

.. autofunction:: tudatpy.dynamics.propagation_setup.integrator.step_size_control_custom_blockwise_scalar_tolerance

.. autofunction:: tudatpy.dynamics.propagation_setup.integrator.step_size_control_custom_blockwise_matrix_tolerance

.. autofunction:: tudatpy.dynamics.propagation_setup.integrator.standard_cartesian_state_element_blocks





Enumerations
------------
.. currentmodule:: tudatpy.dynamics.propagation_setup.integrator

.. autosummary::

   AvailableIntegrators

   CoefficientSets

   OrderToIntegrate

   ExtrapolationMethodStepSequences

   MinimumIntegrationTimeStepHandling



.. autoclass:: tudatpy.dynamics.propagation_setup.integrator.AvailableIntegrators
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.integrator.CoefficientSets
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.integrator.OrderToIntegrate
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.integrator.ExtrapolationMethodStepSequences
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.integrator.MinimumIntegrationTimeStepHandling
   :members:




Classes
-------
.. currentmodule:: tudatpy.dynamics.propagation_setup.integrator

.. autosummary::

   IntegratorStepSizeControlSettings

   IntegratorStepSizeValidationSettings

   IntegratorSettings

   RungeKuttaFixedStepSizeSettings

   BulirschStoerIntegratorSettings

   AdamsBashforthMoultonSettings



.. autoclass:: tudatpy.dynamics.propagation_setup.integrator.IntegratorStepSizeControlSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.integrator.IntegratorStepSizeValidationSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.integrator.IntegratorSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.integrator.RungeKuttaFixedStepSizeSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.integrator.BulirschStoerIntegratorSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.integrator.AdamsBashforthMoultonSettings
   :members:




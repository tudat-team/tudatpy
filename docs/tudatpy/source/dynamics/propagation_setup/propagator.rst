.. _propagator:

``propagator``
==============

This module contains a set of functions for setting up the
propagator settings and associated termination conditions in a propagation.

The propagator settings object(s) (derived classed of :class:`~tudatpy.dynamics.propagation_setup.propagator.PropagatorSettings`)
defined through this module are the container through
which all information on the numerical propagation is communicated to the propagation of dynamics
(through the :func:`~tudatpy.dynamics.simulator.create_dynamics_simulator` function), the propagation of variational equations
(through the :func:`~tudatpy.dynamics.simulator.create_variational_equations_solver` function) or an estimation (through the
:class:`~~tudatpy.estimation.estimation_analysis.Estimator` class constructor).

The propagator settings defined here can be used for propagating single- and multi-arc dynamics, translational dynamics,
rotational dynamics, multi-type dynamics, etc. with a large diversity of termination settings options,
dependent variables to save and many more options. The full options and structure are described on a `dedicated page
<https://docs.tudat.space/en/latest/user-guide/state-propagation/propagation-setup.html>`_

In addition to settings for the full propagation, this module also contains function to create settings for the termination
of a propagation (stored in instances of :class:`~tudatpy.dynamics.propagation_setup.propagator.PropagationTerminationSettings
and its derived class). These objects are passed to the propagator settings function in this module to define the
propagation termination.


.. References
.. ----------
.. .. [1] Vittaldev, V., Mooij, E., & Naeije, M. C. (2012). Unified State Model theory and
..        application in Astrodynamics. Celestial Mechanics and Dynamical Astronomy, 112(3), 253-282.
.. .. [2] Wakker, K. F. (2015). Fundamentals of astrodynamics.
.. .. [3] Hintz, G. R. (2008). Survey of orbit element sets. Journal of guidance, control, and dynamics, 31(3), 785-790.
.. .. [4] Vallado, D. A. (2001). Fundamentals of astrodynamics and applications (Vol. 12). Springer Science & Business Media.






Functions
---------
.. currentmodule:: tudatpy.dynamics.propagation_setup.propagator

.. autosummary::

   translational

   rotational

   mass

   custom_state

   multitype

   multi_arc

   hybrid_arc

   time_termination

   cpu_time_termination

   dependent_variable_termination

   custom_termination

   hybrid_termination

   non_sequential_termination

   add_dependent_variable_settings



.. autofunction:: tudatpy.dynamics.propagation_setup.propagator.translational

.. autofunction:: tudatpy.dynamics.propagation_setup.propagator.rotational

.. autofunction:: tudatpy.dynamics.propagation_setup.propagator.mass

.. autofunction:: tudatpy.dynamics.propagation_setup.propagator.custom_state

.. autofunction:: tudatpy.dynamics.propagation_setup.propagator.multitype

.. autofunction:: tudatpy.dynamics.propagation_setup.propagator.multi_arc

.. autofunction:: tudatpy.dynamics.propagation_setup.propagator.hybrid_arc

.. autofunction:: tudatpy.dynamics.propagation_setup.propagator.time_termination

.. autofunction:: tudatpy.dynamics.propagation_setup.propagator.cpu_time_termination

.. autofunction:: tudatpy.dynamics.propagation_setup.propagator.dependent_variable_termination

.. autofunction:: tudatpy.dynamics.propagation_setup.propagator.custom_termination

.. autofunction:: tudatpy.dynamics.propagation_setup.propagator.hybrid_termination

.. autofunction:: tudatpy.dynamics.propagation_setup.propagator.non_sequential_termination

.. autofunction:: tudatpy.dynamics.propagation_setup.propagator.add_dependent_variable_settings




Enumerations
------------
.. currentmodule:: tudatpy.dynamics.propagation_setup.propagator

.. autosummary::

   TranslationalPropagatorType

   RotationalPropagatorType

   StateType

   PropagationTerminationTypes



.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.TranslationalPropagatorType
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.RotationalPropagatorType
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.StateType
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.PropagationTerminationTypes
   :members:




Classes
-------
.. currentmodule:: tudatpy.dynamics.propagation_setup.propagator

.. autosummary::

   PropagatorSettings

   MultiArcPropagatorSettings

   HybridArcPropagatorSettings

   SingleArcPropagatorSettings

   TranslationalStatePropagatorSettings

   RotationalStatePropagatorSettings

   MultiTypePropagatorSettings

   PropagationTerminationSettings

   PropagationDependentVariableTerminationSettings

   PropagationTimeTerminationSettings

   PropagationCPUTimeTerminationSettings

   PropagationCustomTerminationSettings

   PropagationHybridTerminationSettings

   PropagationPrintSettings

   PropagatorProcessingSettings

   SingleArcPropagatorProcessingSettings

   MultiArcPropagatorProcessingSettings

   HybridArcPropagatorProcessingSettings



.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.PropagatorSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.MultiArcPropagatorSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.HybridArcPropagatorSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.SingleArcPropagatorSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.TranslationalStatePropagatorSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.RotationalStatePropagatorSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.MultiTypePropagatorSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.PropagationTerminationSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.PropagationDependentVariableTerminationSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.PropagationTimeTerminationSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.PropagationCPUTimeTerminationSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.PropagationCustomTerminationSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.PropagationHybridTerminationSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.PropagationPrintSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.PropagatorProcessingSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.SingleArcPropagatorProcessingSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.MultiArcPropagatorProcessingSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.propagator.HybridArcPropagatorProcessingSettings
   :members:




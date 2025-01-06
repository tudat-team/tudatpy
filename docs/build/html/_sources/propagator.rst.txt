``propagator``
==============
This module provides the functionality for creating propagator settings.








References
----------
.. [1] Vittaldev, V., Mooij, E., & Naeije, M. C. (2012). Unified State Model theory and
       application in Astrodynamics. Celestial Mechanics and Dynamical Astronomy, 112(3), 253-282.
.. [2] Wakker, K. F. (2015). Fundamentals of astrodynamics.
.. [3] Hintz, G. R. (2008). Survey of orbit element sets. Journal of guidance, control, and dynamics, 31(3), 785-790.
.. [4] Vallado, D. A. (2001). Fundamentals of astrodynamics and applications (Vol. 12). Springer Science & Business Media.






Functions
---------
.. currentmodule:: tudatpy.numerical_simulation.propagation_setup.propagator

.. autosummary::

   translational

   rotational

   mass

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



.. autofunction:: tudatpy.numerical_simulation.propagation_setup.propagator.translational

.. autofunction:: tudatpy.numerical_simulation.propagation_setup.propagator.rotational

.. autofunction:: tudatpy.numerical_simulation.propagation_setup.propagator.mass

.. autofunction:: tudatpy.numerical_simulation.propagation_setup.propagator.multitype

.. autofunction:: tudatpy.numerical_simulation.propagation_setup.propagator.multi_arc

.. autofunction:: tudatpy.numerical_simulation.propagation_setup.propagator.hybrid_arc

.. autofunction:: tudatpy.numerical_simulation.propagation_setup.propagator.time_termination

.. autofunction:: tudatpy.numerical_simulation.propagation_setup.propagator.cpu_time_termination

.. autofunction:: tudatpy.numerical_simulation.propagation_setup.propagator.dependent_variable_termination

.. autofunction:: tudatpy.numerical_simulation.propagation_setup.propagator.custom_termination

.. autofunction:: tudatpy.numerical_simulation.propagation_setup.propagator.hybrid_termination

.. autofunction:: tudatpy.numerical_simulation.propagation_setup.propagator.non_sequential_termination

.. autofunction:: tudatpy.numerical_simulation.propagation_setup.propagator.add_dependent_variable_settings




Enumerations
------------
.. currentmodule:: tudatpy.numerical_simulation.propagation_setup.propagator

.. autosummary::

   TranslationalPropagatorType

   RotationalPropagatorType

   StateType

   PropagationTerminationTypes



.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.TranslationalPropagatorType
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.RotationalPropagatorType
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.StateType
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.PropagationTerminationTypes
   :members:




Classes
-------
.. currentmodule:: tudatpy.numerical_simulation.propagation_setup.propagator

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



.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.PropagatorSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.MultiArcPropagatorSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.HybridArcPropagatorSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.SingleArcPropagatorSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.TranslationalStatePropagatorSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.RotationalStatePropagatorSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.MultiTypePropagatorSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.PropagationTerminationSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.PropagationDependentVariableTerminationSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.PropagationTimeTerminationSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.PropagationCPUTimeTerminationSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.PropagationCustomTerminationSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.PropagationHybridTerminationSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.PropagationPrintSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.PropagatorProcessingSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.SingleArcPropagatorProcessingSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.MultiArcPropagatorProcessingSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation_setup.propagator.HybridArcPropagatorProcessingSettings
   :members:




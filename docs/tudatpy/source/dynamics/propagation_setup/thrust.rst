.. _thrust:

``thrust``
==========
This module provides the functionality for creating thrust magnitude settings. These models
are used in the calculation of thrust accelerations. The various details of thrust modelling and the interaction
of the functionality in this model with the rest of Tudat is provided on the user guide page for `thrust models <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/thrust_models.html>`_



Functions
---------
.. currentmodule:: tudatpy.dynamics.propagation_setup.thrust

.. autosummary::

   constant_thrust_magnitude

   custom_thrust_magnitude

   custom_thrust_magnitude_fixed_isp

   custom_thrust_acceleration_magnitude

   custom_thrust_acceleration_magnitude_fixed_isp



.. autofunction:: tudatpy.dynamics.propagation_setup.thrust.constant_thrust_magnitude

.. autofunction:: tudatpy.dynamics.propagation_setup.thrust.custom_thrust_magnitude

.. autofunction:: tudatpy.dynamics.propagation_setup.thrust.custom_thrust_magnitude_fixed_isp

.. autofunction:: tudatpy.dynamics.propagation_setup.thrust.custom_thrust_acceleration_magnitude

.. autofunction:: tudatpy.dynamics.propagation_setup.thrust.custom_thrust_acceleration_magnitude_fixed_isp




Enumerations
------------
.. currentmodule:: tudatpy.dynamics.propagation_setup.thrust

.. autosummary::

   ThrustFrames

   ThrustMagnitudeTypes



.. autoclass:: tudatpy.dynamics.propagation_setup.thrust.ThrustFrames
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.thrust.ThrustMagnitudeTypes
   :members:




Classes
-------
.. currentmodule:: tudatpy.dynamics.propagation_setup.thrust

.. autosummary::

   ThrustMagnitudeSettings

   ConstantThrustMagnitudeSettings

   CustomThrustMagnitudeSettings



.. autoclass:: tudatpy.dynamics.propagation_setup.thrust.ThrustMagnitudeSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.thrust.ConstantThrustMagnitudeSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.thrust.CustomThrustMagnitudeSettings
   :members:




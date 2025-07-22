.. _acceleration:

``acceleration``
================

Here, you will find a list of all acceleration models available in Tudat. The functions here all create settings
objects for accelerations that are provided as input to the :func:`~tudatpy.dynamics.propagation_setup.create_acceleration_models`
function. For more details on how these accelerations interface with the rest of the Tudat propagation framework,
see our user guide page on `translational dynamics <https://docs.tudat.space/en/latest/user-guide/state-propagation/propagation-setup/translational.html>`_

Mathematical model descriptions are provided either directly here per acceleration, or on a dedicated page of the Tudat user guide for
 * `Gravitational accelerations <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/third_body_acceleration.html>`_
 * `Thrust accelerations <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/thrust_models.html>`_
 * `Radiation pressure accelerations <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/radiation_pressure_acceleration.html>`_.

The functions in this submodule create objects of type :class:`~AccelerationSettings` or (in case more information than only
the type of acceleration is needed to create the acceleration), one of its derived classes.

Functions
---------
.. currentmodule:: tudatpy.dynamics.propagation_setup.acceleration

.. autosummary::

   point_mass_gravity

   spherical_harmonic_gravity

   mutual_spherical_harmonic_gravity

   aerodynamic

   radiation_pressure

   polyhedron_gravity

   ring_gravity

   relativistic_correction

   einstein_infeld_hofmann

   empirical

   custom_acceleration

   direct_tidal_dissipation_acceleration

   quasi_impulsive_shots_acceleration

   thrust_from_engines

   thrust_from_engine

   thrust_from_all_engines

   yarkovsky



.. autofunction:: tudatpy.dynamics.propagation_setup.acceleration.point_mass_gravity

.. autofunction:: tudatpy.dynamics.propagation_setup.acceleration.spherical_harmonic_gravity

.. autofunction:: tudatpy.dynamics.propagation_setup.acceleration.mutual_spherical_harmonic_gravity

.. autofunction:: tudatpy.dynamics.propagation_setup.acceleration.aerodynamic

.. autofunction:: tudatpy.dynamics.propagation_setup.acceleration.radiation_pressure

.. autofunction:: tudatpy.dynamics.propagation_setup.acceleration.polyhedron_gravity

.. autofunction:: tudatpy.dynamics.propagation_setup.acceleration.ring_gravity

.. autofunction:: tudatpy.dynamics.propagation_setup.acceleration.einstein_infeld_hofmann

.. autofunction:: tudatpy.dynamics.propagation_setup.acceleration.relativistic_correction

.. autofunction:: tudatpy.dynamics.propagation_setup.acceleration.empirical

.. autofunction:: tudatpy.dynamics.propagation_setup.acceleration.custom_acceleration

.. autofunction:: tudatpy.dynamics.propagation_setup.acceleration.direct_tidal_dissipation_acceleration

.. autofunction:: tudatpy.dynamics.propagation_setup.acceleration.quasi_impulsive_shots_acceleration

.. autofunction:: tudatpy.dynamics.propagation_setup.acceleration.thrust_from_engines

.. autofunction:: tudatpy.dynamics.propagation_setup.acceleration.thrust_from_engine

.. autofunction:: tudatpy.dynamics.propagation_setup.acceleration.thrust_from_all_engines

.. autofunction:: tudatpy.dynamics.propagation_setup.acceleration.yarkovsky




Enumerations
------------
.. currentmodule:: tudatpy.dynamics.propagation_setup.acceleration

.. autosummary::

   AvailableAcceleration



.. autoclass:: tudatpy.dynamics.propagation_setup.acceleration.AvailableAcceleration
   :members:




Classes
-------
.. currentmodule:: tudatpy.dynamics.propagation_setup.acceleration

.. autosummary::

   AccelerationSettings

   SphericalHarmonicAccelerationSettings

   MutualSphericalHarmonicAccelerationSettings

   RelativisticAccelerationCorrectionSettings

   EmpiricalAccelerationSettings

   CustomAccelerationSettings

   DirectTidalDissipationAccelerationSettings

   ThrustAccelerationSettings

   MomentumWheelDesaturationAccelerationSettings



.. autoclass:: tudatpy.dynamics.propagation_setup.acceleration.AccelerationSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.acceleration.SphericalHarmonicAccelerationSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.acceleration.MutualSphericalHarmonicAccelerationSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.acceleration.RelativisticAccelerationCorrectionSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.acceleration.EmpiricalAccelerationSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.acceleration.CustomAccelerationSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.acceleration.DirectTidalDissipationAccelerationSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.acceleration.ThrustAccelerationSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.acceleration.MomentumWheelDesaturationAccelerationSettings
   :members:
``ephemeris``
=============
This module contains a set of factory functions for setting up the
ephemeris models of celestial bodies in an environment. Below a short
overview of aspects of some of the ephemeris models in order to aid in
properly selecting an choosing a model.

**Spice-based models** For many typical applications, natural body ephemerides
will be calculated from `Spice kernels <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/default_env_models.html#spice-in-tudat>`_.
In some cases, a user may find that the default Spice kernels are insufficient
for their purposes, due to one of two reasons:

* The body for which the state is required *is* in the ephemeris Spice kernel, but the time at which the state is needed lies outside of the bounds for which the Spice kernel has data
* The body for which the state is required *is not* in the ephemeris Spice kernel

In both cases, a user should load additional Spice kernels. This can be done using the :func:`~tudatpy.interface.spice.load_kernel`. Spice kernels for many bodies may be found in a number of places.
The 'goto' place for Spice kernels for ephemerides is the NAIF website (developers of Spice), which you can find
`here <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/>`_.

**Use of scaled models** For a sensitivity analysis (among others) it may be useful to modify the ephemeris of a body, for instance
to emulate the influence of a 1 km offset in the state provided by the nominal ephemeris. Unlike most other environment models,
this cannot be achieved (at least not for most types of ephemerides) by modifying a single defining parameter of the model.
Instead, we provide the functions
:func:`~tudatpy.numerical_simulation.environment_setup.ephemeris.scaled_by_vector` and
:func:`~tudatpy.numerical_simulation.environment_setup.ephemeris.scaled_by_vector_function`,
which take nominal ephemeris settings, and add a user-defined variation (constant or time-varying; absolute or relative) to the
inertial Cartesian state elements produced by the ephemeris.

**Using the ephemeris outside the propagation** In various cases, the ephemeris object is useful to use independently of the propagation. Details can be found in the API entry for :class:`~tudatpy.numerical_simulation.environment.Ephemeris`, but we provide a short example here as well.

.. code-block:: python

  bodies = .... // Create system of bodies
  earth_ephemeris = bodies.get('Earth').ephemeris
  earth_state_at_epoch = earth_ephemeris.cartesian_state( epoch )

where the ``epoch`` input is (as always in Tudat) the time in seconds since J2000. The ``earth_state_at_epoch`` is always in a frame with inertial orientation. The specific orientation and origin can be access from the :attr:`~tudatpy.numerical_simulation.environment.Ephemeris.frame_orientation` and :attr:`~tudatpy.numerical_simulation.environment.Ephemeris.frame_origin` attributes.

* jpl_horizons (:func:`~tudatpy.numerical_simulation.environment_setup.ephemeris.horizons_wrapper.jpl_horizons`)










Functions
---------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.ephemeris

.. autosummary::

   direct_spice

   interpolated_spice

   approximate_jpl_model

   constant

   custom_ephemeris

   keplerian

   keplerian_from_spice

   scaled_by_constant

   scaled_by_vector

   scaled_by_vector_function

   tabulated

   tabulated_from_existing

   horizons_wrapper.jpl_horizons


.. autofunction:: tudatpy.numerical_simulation.environment_setup.ephemeris.direct_spice

.. autofunction:: tudatpy.numerical_simulation.environment_setup.ephemeris.interpolated_spice

.. autofunction:: tudatpy.numerical_simulation.environment_setup.ephemeris.approximate_jpl_model

.. autofunction:: tudatpy.numerical_simulation.environment_setup.ephemeris.constant

.. autofunction:: tudatpy.numerical_simulation.environment_setup.ephemeris.custom_ephemeris

.. autofunction:: tudatpy.numerical_simulation.environment_setup.ephemeris.keplerian

.. autofunction:: tudatpy.numerical_simulation.environment_setup.ephemeris.keplerian_from_spice

.. autofunction:: tudatpy.numerical_simulation.environment_setup.ephemeris.scaled_by_constant

.. autofunction:: tudatpy.numerical_simulation.environment_setup.ephemeris.scaled_by_vector

.. autofunction:: tudatpy.numerical_simulation.environment_setup.ephemeris.scaled_by_vector_function

.. autofunction:: tudatpy.numerical_simulation.environment_setup.ephemeris.tabulated

.. autofunction:: tudatpy.numerical_simulation.environment_setup.ephemeris.tabulated_from_existing

.. autofunction:: tudatpy.numerical_simulation.environment_setup.ephemeris.horizons_wrapper.jpl_horizons




Classes
-------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.ephemeris

.. autosummary::

   EphemerisSettings

   ScaledEphemerisSettings

   DirectSpiceEphemerisSettings

   InterpolatedSpiceEphemerisSettings

   ApproximateJplEphemerisSettings

   ConstantEphemerisSettings

   CustomEphemerisSettings

   TabulatedEphemerisSettings



.. autoclass:: tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.ephemeris.ScaledEphemerisSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.ephemeris.DirectSpiceEphemerisSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.ephemeris.InterpolatedSpiceEphemerisSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.ephemeris.ApproximateJplEphemerisSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.ephemeris.ConstantEphemerisSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.ephemeris.CustomEphemerisSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.ephemeris.TabulatedEphemerisSettings
   :members:




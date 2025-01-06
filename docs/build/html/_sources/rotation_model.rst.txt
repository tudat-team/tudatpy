``rotation_model``
==================
This module contains a set of factory functions for setting up the
rotation models of celestial bodies in an environment. Below a short
overview of aspects of some of the rotation models in order to aid in
properly selecting an choosing a model.

Tudat has a broad range of rotation models available. In principle, these models can be assigned to both celestial bodies and natural bodies. 
However, a subset of these models is typically only applied to natural *or* artificial bodies. Rotation models have a wide range of,
sometimes indirect, influences on the dynamics

* A spherical harmonic acceleration exerted by a central body is first evaluated in a body-fixed frame, and the transformed to an inertial frame. Consequently, the central body's rotation has a fundamental influence on the exerted spherical harmonic acceleration
* A `thrust acceleration <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/thrust_models.html#thrust-models>`_ in Tudat is calculated from two models: (1) an engine model, which defined the body-fixed direction of the thrust, and the magnitude of the thrust (2) the orientation of the body in space, defined by its rotation model
* For a non-spherical central body shape models, the current orientation of this central body has an indirect influence on the altitude at which a vehicle with a given *inertial* state is located

**Rotation and thrust** Two rotation models, which are typically used for vehicles under `thrust <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/thrust_models.html#thrust-models>`_, and/or vehicles undergoing `aerodynamic forces <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/aerodynamics.html#aerodynamic-models>`_, are the following:

* The rotation model :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.aerodynamic_angle_based`, which calculates the body's rotation based on the angle of attack, sideslip angle and bank angle. Note that these angles are definend w.r.t. the relative wind. This model is typical when using, for instance, a re-entry simulation. It imposes these three angles, and calculates the body orientation by combination with the latitude, longitude, heading angle, flight path angles. There is a related model, :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.zero_pitch_moment_aerodynamic_angle_based`, that uses the same setup, but does not impose the angle of attack, but caculates by imposing aerodynamic pitch trim (zero pitch moment).
* The rotation model :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.custom_inertial_direction_based`, which is typical when calculating dynamics of a vehicle under thrust. It is based on linking a body-fixed  direction (now limited to the body-fixed x-axis) to an arbitrary inertial direction. This allows the thrust (assuming that this is aligned with this same body-fixed direction) to be guided in an inertial direction determined by a user-defined model. 

**Relation to gravity field** When modifying the rotation model settings, the name of the body-fixed frame may also be changed (as is the case for, for instance, the :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.gcrs_to_itrs`, where the body-fixed frame has the name "ITRS").
One consequence of this is that you may get an error from the spherical harmonic gravity field, which can no longer find the frame to which it is associated. This can be resolved by (for instance) associating the gravity field to the new frame. For the above example, this would be done by the following:

.. code-block:: python
                
    body_settings.get( "Earth" ).gravity_field_settings.associated_reference_frame = "ITRS"
    
**High-accuracy Earth rotation model** The :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.gcrs_to_itrs` creates a high accuracy rotation model, following the IERS 2010 Conventions. This includes small variations that are not predicted by models, but are instead measured by geodetic techniques and published as tabulated data by the IERS. If so desired, the exact files used for these corrections may be adapted by the user (see :func:`~tudatpy.astro.earth_orientation.EarthOrientationAnglesCalculator`), which includes specific settings for daily variations in earth rotation angle, which influences the UTC - UT1 time conversion. 

**Using the rotation model outside the propagation** In various cases, the rotation model object is useful to use independently of the propagation. Details can be found in the API entry for :class:`~tudatpy.numerical_simulation.environment.RotationalEphemeris`, but we provide a short example here as well.

.. code-block:: python

    bodies = .... // Create system of bodies
    earth_rotation_model = bodies.get('Earth').rotation_model
    earth_rotation_at_epoch = earth_rotation_model.body_fixed_to_inertial_rotation( epoch )

where the ``epoch`` input is (as always in Tudat) the time in seconds since J2000. The specific rotation model provides the orientation from the :attr:`~tudatpy.numerical_simulation.environment.RotationalEphemeris.inertial_frame_name` to the :attr:`~tudatpy.numerical_simulation.environment.RotationalEphemeris.body_fixed_frame_name` frames. In the above example, the rotation matrix from the body-fixed to the inertial frame is extracted. Other functions are available in the :class:`~tudatpy.numerical_simulation.environment.RotationalEphemeris` to extract the inverse rotation, its time-derivative, and the angular velocity vector of the body-fixed frame. Finally, note that the :func:`~tudatpy.numerical_simulation.environment.transform_to_inertial_orientation`, which uses the rotation model to rotation a body-fixed to an inertial state, may be useful in this context for some applications.












Functions
---------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.rotation_model

.. autosummary::

   simple

   simple_from_spice

   synchronous

   spice

   gcrs_to_itrs

   constant_rotation_model

   aerodynamic_angle_based

   zero_pitch_moment_aerodynamic_angle_based

   custom_inertial_direction_based

   orbital_state_direction_based

   mars_high_accuracy



.. autofunction:: tudatpy.numerical_simulation.environment_setup.rotation_model.simple

.. autofunction:: tudatpy.numerical_simulation.environment_setup.rotation_model.simple_from_spice

.. autofunction:: tudatpy.numerical_simulation.environment_setup.rotation_model.synchronous

.. autofunction:: tudatpy.numerical_simulation.environment_setup.rotation_model.spice

.. autofunction:: tudatpy.numerical_simulation.environment_setup.rotation_model.gcrs_to_itrs

.. autofunction:: tudatpy.numerical_simulation.environment_setup.rotation_model.constant_rotation_model

.. autofunction:: tudatpy.numerical_simulation.environment_setup.rotation_model.aerodynamic_angle_based

.. autofunction:: tudatpy.numerical_simulation.environment_setup.rotation_model.zero_pitch_moment_aerodynamic_angle_based

.. autofunction:: tudatpy.numerical_simulation.environment_setup.rotation_model.custom_inertial_direction_based

.. autofunction:: tudatpy.numerical_simulation.environment_setup.rotation_model.orbital_state_direction_based

.. autofunction:: tudatpy.numerical_simulation.environment_setup.rotation_model.mars_high_accuracy




Enumerations
------------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.rotation_model

.. autosummary::

   RotationModelType

   IAUConventions



.. autoclass:: tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelType
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.rotation_model.IAUConventions
   :members:




Classes
-------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.rotation_model

.. autosummary::

   RotationModelSettings



.. autoclass:: tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings
   :members:




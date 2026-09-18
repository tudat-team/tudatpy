.. _gravity_field:

``gravity_field``
=================
This module contains a set of factory functions for setting up the
gravitational potential models of celestial bodies in an environment.

The main interfaces with Tudat is the :attr:`~tudatpy.dynamics.environment_setup.BodySettings.gravity_field_settings`
attribute  (of type :class:`~tudatpy.dynamics.environment_setup.gravity_field.GravityFieldSettings`) of the body settings, which defines settings for the gravity field of a body.
**The functions in this submodule are used to create these settings objects.** When creating a body (typically using the
:func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies` function), an object of type
:class:`~tudatpy.dynamics.environment.GravityFieldModel` (or a derived class) is created
and added to the associated :class:`~tudatpy.dynamics.environment.Body` object based on the settings object, which can
be retrieved using the :attr:`~tudatpy.dynamics.environment.Body.gravity_field_model` attribute.

The gravity field model is used in numerous aspects of a propagation and estimation in Tudat. Most prominently, it
is used to compute gravitational accelerations (see :ref:`acceleration`) and torques (see :ref:`torque`).

The following code block gives an overview of the steps to define, create, and extract a gravity field model, for the specific example of a point-mass model with :math:`\mu=3.986004418\cdot 10^{14}` m\ :sup:`3`/s\ :sup:`2`.

.. code-block:: python

  from tudatpy.dynamics import environment_setup

  # Create body settings
  body_settings =  environment_setup.get_default_body_settings( ... ) # Typical way to instantiate body settings

  # Modify gravity field model settings (base class type GravityFieldSettings)
  body_settings.get( 'Earth' ).gravity_field_settings = environment_setup.gravity_field.central(
      gravitational_parameter = 3.986004418E14 )

  # Create bodies
  bodies = environment_setup.create_system_of_bodies(body_settings)

  # Extract gravity field model (base class type GravityFieldModel) from Earth
  earth_gravity_field_model = bodies.get( 'Earth' ).gravity_field_model

Below a short
overview of aspects of some of the gravity field models in order to aid in
properly selecting an choosing a model.
Unlike most other environment model options in Tudat, there are multiple options for creating either a spherical harmonic gravity field, and a point mass gravity field:

* Point-mass gravity field: defining the gravitational parameter manually (:func:`~tudatpy.dynamics.environment_setup.gravity_field.central`), extracting it from Spice (:func:`~tudatpy.dynamics.environment_setup.gravity_field.central_spice`) or from the Small-Body Database (SBDB; :func:`~tudatpy.dynamics.environment_setup.gravity_field.sbdb_wrapper.central_sbdb`)
* Spherical harmonic gravity field: defining all the settings manually (:func:`~tudatpy.dynamics.environment_setup.gravity_field.spherical_harmonic`), loading a pre-defined model for a solar system body (:func:`~tudatpy.dynamics.environment_setup.gravity_field.from_file_spherical_harmonic`) or calculating the spherical harmonic coefficients (up to a given degree) based on an ellipsoidal homogeneous mass distribution (:func:`~tudatpy.dynamics.environment_setup.gravity_field.sh_triaxial_ellipsoid_from_density` and :func:`~tudatpy.dynamics.environment_setup.gravity_field.sh_triaxial_ellipsoid_from_gravitational_parameter`)

In Tudat, the gravity field itself is not directly responsible for owning the mass, center of mass or inertia tensor used by a simulation. This is handled by 'rigid body properties' (see :ref:`rigid_body`). Unless explicit rigid-body settings already exist, compatible gravity-derived settings are created automatically in the following manner:

* Point-mass gravity field: mass computed from gravitational parameter, no inertia tensor, and center of mass at origin of body-fixed frame
* Spherical harmonic gravity field: mass computed from gravitational parameter, center of mass computed from degree-one coefficients, and inertia computed from complete degree-two coefficients plus the scaled mean moment configured with :func:`~tudatpy.dynamics.environment_setup.rigid_body.from_gravity_field`. A gravity field remains valid when no inertia tensor can be defined.
* Polyhedron gravity field: mass computed from gravitational parameter and inertia tensor computed from the homogeneous mass distribution represented by the polyhedron.
* Ring gravity field model: mass computed from gravitational parameter, center of mass at the body-fixed origin, and no inertia tensor.

Gravity models retain a link to the associated rigid-body properties so that changes to the gravitational parameter or relevant spherical-harmonic coefficients immediately refresh gravity-derived mass, center of mass, and inertia. Explicit non-gravity rigid-body settings retain precedence and are not overwritten by these callbacks.



.. References
.. ----------
.. .. [1] Balmino, G. (1994). Gravitational potential harmonics from the shape of an homogeneous body. Celestial
..       Mechanics and Dynamical Astronomy, 60(3), 331-364.
.. .. [2] Werner, R. A., and Scheeres, D. J. (1997). Exterior Gravitation of a Polyhedron Derived and Compared With
..       Harmonic and Mascon Gravitation Representations of Asteroid 4769 Castalia. Celestial Mechanics and Dynamical
..       Astronomy, 65, 313-344.






Functions
---------
.. currentmodule:: tudatpy.dynamics.environment_setup.gravity_field

.. autosummary::

   central

   central_spice

   spherical_harmonic

   sh_triaxial_ellipsoid_from_density

   sh_triaxial_ellipsoid_from_gravitational_parameter

   from_file_spherical_harmonic

   predefined_spherical_harmonic

   polyhedron_from_mu

   polyhedron_from_density

   ring_model

   sbdb_wrapper.central_sbdb

   sbdb_wrapper.central_sbdb_density


.. autofunction:: tudatpy.dynamics.environment_setup.gravity_field.central

.. autofunction:: tudatpy.dynamics.environment_setup.gravity_field.central_spice

.. autofunction:: tudatpy.dynamics.environment_setup.gravity_field.spherical_harmonic

.. autofunction:: tudatpy.dynamics.environment_setup.gravity_field.sh_triaxial_ellipsoid_from_density

.. autofunction:: tudatpy.dynamics.environment_setup.gravity_field.sh_triaxial_ellipsoid_from_gravitational_parameter

.. autofunction:: tudatpy.dynamics.environment_setup.gravity_field.from_file_spherical_harmonic

.. autofunction:: tudatpy.dynamics.environment_setup.gravity_field.predefined_spherical_harmonic

.. autofunction:: tudatpy.dynamics.environment_setup.gravity_field.polyhedron_from_mu

.. autofunction:: tudatpy.dynamics.environment_setup.gravity_field.polyhedron_from_density

.. autofunction:: tudatpy.dynamics.environment_setup.gravity_field.ring_model

.. autofunction:: tudatpy.dynamics.environment_setup.gravity_field.sbdb_wrapper.central_sbdb

.. autofunction:: tudatpy.dynamics.environment_setup.gravity_field.sbdb_wrapper.central_sbdb_density


Enumerations
------------
.. currentmodule:: tudatpy.dynamics.environment_setup.gravity_field

.. autosummary::

   GravityFieldType

   PredefinedSphericalHarmonicsModel



.. autoclass:: tudatpy.dynamics.environment_setup.gravity_field.GravityFieldType
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.gravity_field.PredefinedSphericalHarmonicsModel
   :members:




Classes
-------
.. currentmodule:: tudatpy.dynamics.environment_setup.gravity_field

.. autosummary::

   GravityFieldSettings

   CentralGravityFieldSettings

   SphericalHarmonicsGravityFieldSettings

   PolyhedronGravityFieldSettings



.. autoclass:: tudatpy.dynamics.environment_setup.gravity_field.GravityFieldSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.gravity_field.CentralGravityFieldSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.gravity_field.PolyhedronGravityFieldSettings
   :members:

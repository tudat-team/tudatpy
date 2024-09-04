``gravity_field``
=================
This module contains a set of factory functions for setting up the
gravitational potential models of celestial bodies in an environment. Below a short
overview of aspects of some of the gravity field models in order to aid in
properly selecting an choosing a model.
Unlike most other environment model options in Tudat, there are multiple options for creating either a spherical harmonic gravity field, and a point mass gravity field:

* Point-mass gravity field: defining the gravitational parameter manually (:func:`~tudatpy.numerical_simulation.environment_setup.gravity_field.central`) or requiring the gravitational parameter to be extracted from Spice (:func:`~tudatpy.numerical_simulation.environment_setup.gravity_field.central_spice`).
* Spherical harmonic gravity field: defining all the settings manually (:func:`~tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic`), loading a pre-defined model for a solar system body (:func:`~tudatpy.numerical_simulation.environment_setup.gravity_field.from_file_spherical_harmonic`) or calculating the spherical harmonic coefficients (up to a given degree) based on an ellipsoidal homogeneous mass distribution (:func:`~tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic_triaxial_body`)

Rigid body properties will always be created automatically when a body is endowed with a gravity field, as described below:

* Point-mass gravity field: mass computed from gravitational parameter; zero inertia tensor, and center of mass at origin of body-fixed frame
* Spherical harmonic gravity field: mass computed from gravitational parameter, center of mass computed from degree 1 gravity field coefficients, inertia tensor as described in :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic`
* Polyhedron gravity field: mass computed from gravitational parameter, center of mass and inertia tensor computed from homogeneous mas distribution inside body

* sbdb_wrapper (:func:`~tudatpy.numerical_simulation.environment_setup.gravity_field.sbdb_wrapper.central_sbdb`)






References
----------
.. [1] Balmino, G. (1994). Gravitational potential harmonics from the shape of an homogeneous body. Celestial
      Mechanics and Dynamical Astronomy, 60(3), 331-364.
.. [2] Werner, R. A., and Scheeres, D. J. (1997). Exterior Gravitation of a Polyhedron Derived and Compared With
      Harmonic and Mascon Gravitation Representations of Asteroid 4769 Castalia. Celestial Mechanics and Dynamical
      Astronomy, 65, 313-344.






Functions
---------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.gravity_field

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

   sbdb_wrapper.central_sbdb


.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field.central

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field.central_spice

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field.sh_triaxial_ellipsoid_from_density

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field.sh_triaxial_ellipsoid_from_gravitational_parameter

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field.from_file_spherical_harmonic

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field.predefined_spherical_harmonic

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field.polyhedron_from_mu

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field.polyhedron_from_density

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field.sbdb_wrapper.central_sbdb


Enumerations
------------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.gravity_field

.. autosummary::

   GravityFieldType

   PredefinedSphericalHarmonicsModel



.. autoclass:: tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldType
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.gravity_field.PredefinedSphericalHarmonicsModel
   :members:




Classes
-------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.gravity_field

.. autosummary::

   GravityFieldSettings

   CentralGravityFieldSettings

   SphericalHarmonicsGravityFieldSettings

   PolyhedronGravityFieldSettings



.. autoclass:: tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.gravity_field.CentralGravityFieldSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.gravity_field.PolyhedronGravityFieldSettings
   :members:




.. _shape_deformation:

``shape_deformation``
=====================
This module contains a set of factory functions for setting up the
shape deformation models of celestial bodies in an environment.

The main interfaces with Tudat is the :attr:`~tudatpy.dynamics.environment_setup.BodySettings.shape_deformation_settings` list
attribute  (with entries of type :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BodyShapeSettings`) of the body settings, which defines settings for the global shape deformation model of a body.
**The functions in this submodule are used to create these settings objects.** When creating a body (typically using the
:func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies` function), a list of objects of type
:class:`~tudatpy.dynamics.environment.BodyShapeDeformationModel` (or a derived class) is created
and added to the associated :class:`~tudatpy.dynamics.environment.Body` object based on the settings object, which can
be retrieved using the :attr:`~tudatpy.dynamics.environment.Body.shape_deformation_model` attribute.

The shape deformation models are used for global shape variations of a celestial body, which are used (for instance)
for the high-accuracy modelling of the positions of ground stations on a body. Station-specific models for position variations
of a station (such as plate motion) are set in the ground stations themselves (see :ref:`ground_station`).


Ocean tidal loading evaluates TT seconds since J2000 and uses corresponding UTC phases.
Standalone calls require TT. The existing environment currently passes its TDB time as TT,
retaining the millisecond TDB-TT approximation until explicit time-scale routing is implemented.

Functions
---------
.. currentmodule:: tudatpy.dynamics.environment_setup.shape_deformation

.. autosummary::

   basic_solid_body_tidal

   degree_two_basic_solid_body_tidal

   iers_2010_solid_body_tidal

   ocean_tidal

   ocean_tidal_from_data

   pole_tidal
.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.basic_solid_body_tidal

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.degree_two_basic_solid_body_tidal

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.iers_2010_solid_body_tidal

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.ocean_tidal

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.ocean_tidal_from_data

.. autoclass:: tudatpy.dynamics.environment_setup.shape_deformation.OceanTideSynthesis
   :members:

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.pole_tidal





Conventional ocean pole-tide loading implements IERS Section 7.1.5 with the
2018 secular pole and IERS 2010 Chapter 1 constants. The official grid covers
latitudes -89.75 through +89.75 degrees. ``require_grid_coverage`` raises outside
the supplied extent; ``clamp_to_nearest_grid_latitude_approximation`` explicitly
uses the nearest covered latitude. No polar-cap interpolation is specified by
the IERS coefficient grid. Longitude is periodic for a complete grid. Cartesian
displacement uses the geocentric spherical basis; public station ENU uses its
actual topocentric frame. Standalone epochs are TT; the existing environment
currently approximates its TDB time argument as TT.

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.ocean_pole_tidal

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.ocean_pole_tidal_from_data

.. autoclass:: tudatpy.dynamics.environment_setup.shape_deformation.OceanPoleTideLatitudePolicy
   :members:

Classes
-------
.. currentmodule:: tudatpy.dynamics.environment_setup.shape_deformation

.. autosummary::

   BodyDeformationSettings

   BasicSolidBodyDeformationSettings



.. autoclass:: tudatpy.dynamics.environment_setup.shape_deformation.BodyDeformationSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.shape_deformation.BasicSolidBodyDeformationSettings
   :members:

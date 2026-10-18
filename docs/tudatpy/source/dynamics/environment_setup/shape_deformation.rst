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


In addition to tidal deformation, models for the loading deformation of the Earth are available: ocean tidal loading
(from BLQ files), ocean pole tide loading, atmospheric tidal loading, and precomputed station loading series (for
instance non-tidal atmospheric, oceanic or hydrological loading). The input data of these models are read with the
functions in :ref:`environment_data`.

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

   ocean_pole_tidal

   ocean_pole_tidal_from_data

   atmospheric_tide_from_grid

   atmospheric_tide_from_stations

   atmospheric_tide_from_grid_file

   atmospheric_tide_from_station_file

   station_loading

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.basic_solid_body_tidal

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.degree_two_basic_solid_body_tidal

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.iers_2010_solid_body_tidal

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.ocean_tidal

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.ocean_tidal_from_data

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.pole_tidal

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.ocean_pole_tidal

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.ocean_pole_tidal_from_data

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.atmospheric_tide_from_grid

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.atmospheric_tide_from_stations

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.atmospheric_tide_from_grid_file

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.atmospheric_tide_from_station_file

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.station_loading



Enumerations
------------
.. currentmodule:: tudatpy.dynamics.environment_setup.shape_deformation

.. autosummary::

   OceanTideSynthesis

.. autoclass:: tudatpy.dynamics.environment_setup.shape_deformation.OceanTideSynthesis
   :members:


Classes
-------
.. currentmodule:: tudatpy.dynamics.environment_setup.shape_deformation

.. autosummary::

   BodyDeformationSettings

   BasicSolidBodyDeformationSettings

   OceanPoleTideBodyDeformationSettings

   AtmosphericTideBodyDeformationSettings

   AtmosphericTideEvaluationSettings

   StationLoadingBodyDeformationSettings



.. autoclass:: tudatpy.dynamics.environment_setup.shape_deformation.BodyDeformationSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.shape_deformation.BasicSolidBodyDeformationSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.shape_deformation.OceanPoleTideBodyDeformationSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.shape_deformation.AtmosphericTideBodyDeformationSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.shape_deformation.AtmosphericTideEvaluationSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.shape_deformation.StationLoadingBodyDeformationSettings
   :members:


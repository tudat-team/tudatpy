.. _ground_station:

``ground_station``
==================
This module contains a set of factory functions for setting up ground
stations and associated models. Note that in Tudat, no distinction is
made between a ground station/lander on Earth or a different body.

The main interfaces with Tudat is the :attr:`~tudatpy.dynamics.environment_setup.BodySettings.ground_station_settings` list
attribute (with entries of type :class:`~tudatpy.dynamics.environment_setup.ground_station.GroundStationSettings`) of the body settings, which defines settings for the ground stations of a body.
The functions in this submodule are used to create these settings objects. When creating a body (typically using the
:func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies` function), a list of objects of type
:class:`~tudatpy.dynamics.environment.GroundStation` (or a derived class) is created
and added to the associated :class:`~tudatpy.dynamics.environment.Body` object based on the settings object, which can
be retrieved using the :attr:`~tudatpy.dynamics.environment.Body.ground_station_list` attribute. The ground station
can also be safely added to the ``Body`` after its creation using the :func:`~tudatpy.dynamics.environment.add_ground_station` function.

A ground station defines a reference point (and other relevant properties)
on a celestial body. Although ground stations are considered part of the 
environment in Tudat (as properties of a ``Body`` object), they do not 
influence the numerical propagation (unless a custom model imposing this 
is implemented by the user).

Properties such as a transmitting frequency of the station is added after its creation (to the :class:`~tudatpy.dynamics.environment.GroundStation` object). When creating the settings for a ground stations (using the functions in this module)
the body may be endowed with a list of :class:`~tudatpy.dynamics.environment_setup.ground_stations.GroundStationMotionSettings`
settings, which define station dependent models to define the station-specific motion of a ground station (such as plate motion).
Models for the deformation of the full body (such as tidal shape varitiations) are to be defined through the :ref:`shape_deformation`
module.

For Earth, we provide several options to create default stations, such as the :func:`~tudatpy.dynamics.environment_setup.dsn_stations` and
:func:`~tudatpy.dynamics.environment_setup.evn_stations`.













Functions
---------
.. currentmodule:: tudatpy.dynamics.environment_setup.ground_station

.. autosummary::

    basic_station

    dsn_station

    dsn_stations

    radio_telescope_stations

    evn_stations

    linear_station_motion

    piecewise_constant_station_motion

    custom_station_motion



.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.basic_station

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.dsn_station

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.dsn_stations

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.radio_telescope_stations

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.evn_stations

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.linear_station_motion

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.piecewise_constant_station_motion

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.custom_station_motion






Classes
-------
.. currentmodule:: tudatpy.dynamics.environment_setup.ground_station

.. autosummary::

   GroundStationSettings

   GroundStationMotionSettings

   LinearGroundStationMotionSettings

   PiecewiseConstantGroundStationMotionSettings

   CustomGroundStationMotionSettings



.. autoclass:: tudatpy.dynamics.environment_setup.ground_station.GroundStationSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.ground_station.GroundStationMotionSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.ground_station.LinearGroundStationMotionSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.ground_station.PiecewiseConstantGroundStationMotionSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.ground_station.CustomGroundStationMotionSettings
   :members:




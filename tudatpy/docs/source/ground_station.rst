``ground_station``
==================
This module contains a set of factory functions for setting up ground
stations and associated models. Note that in Tudat, no distinction is
made between a ground station/lander on Earth or a different body. 
A ground station defines a reference point (and other relevant properties)
on a celestial body. Although ground stations are considered part of the 
environment in Tudat (as properties of a ``Body`` object), they do not 
influence the numerical propagation (unless a custom model imposing this 
is implemented by the user). Ground stations can be defined through the 
``BodySettings`` as any other model. But, as the rest of the environment 
does not depend on them, they can safely be added to a body after it is 
created. The process is similar to the one described for :ref:`decorate_empty_body`. 
Specifically, ground station settings are created, and these are then used 
to create a ground station and add it to the body. The specifics of creating
ground station settings is described 
`in the API documentation <https://py.api.tudat.space/en/latest/ground_stations.html>`_. 
An example is given below:

.. code-block:: python 

  # Create ground station settings 
  ground_station_settings = environment_setup.ground_station.basic_station( 
      "TrackingStation", 
      [station_altitude, delft_latitude, delft_longitude], 
      element_conversion.geodetic_position_type) 

  # Add the ground station to the environment 
  environment_setup.add_ground_station( 
      bodies.get_body("Earth"), 
      ground_station_settings ) 


where a simple ground station is created (with only a name and a position), with its position defined in geodetic elements. The position of a ground station in a body-fixed frame can have two sources of time-variability:

* From `shape deformation models <https://py.api.tudat.space/en/latest/shape_deformation.html>`_ of the body on which it is located
* From a list of :class:`~tudatpy.numerical_simulation.environment_setup.ground_station.GroundStationMotionSettings` objects, which can be assigned to the ground station settings (see e.g. :func:`~tudatpy.numerical_simulation.environment_setup.ground_station.basic_station`). These models define time-variability of individual ground stations, in addition to the global shape deformation.

To automatically create a list of settings for all DSN stations (which are then typically assigned to the ``ground_station_settings`` of Earth), the :func:`~tudatpy.numerical_simulation.environment_setup.ground_station.dsn_station_settings` can be used.












Functions
---------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.ground_station

.. autosummary::

   basic_station

   dsn_stations

   linear_station_motion

   piecewise_constant_station_motion

   custom_station_motion



.. autofunction:: tudatpy.numerical_simulation.environment_setup.ground_station.basic_station

.. autofunction:: tudatpy.numerical_simulation.environment_setup.ground_station.dsn_stations

.. autofunction:: tudatpy.numerical_simulation.environment_setup.ground_station.linear_station_motion

.. autofunction:: tudatpy.numerical_simulation.environment_setup.ground_station.piecewise_constant_station_motion

.. autofunction:: tudatpy.numerical_simulation.environment_setup.ground_station.custom_station_motion






Classes
-------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.ground_station

.. autosummary::

   GroundStationSettings

   GroundStationMotionSettings

   LinearGroundStationMotionSettings

   PiecewiseConstantGroundStationMotionSettings

   CustomGroundStationMotionSettings



.. autoclass:: tudatpy.numerical_simulation.environment_setup.ground_station.GroundStationSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.ground_station.GroundStationMotionSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.ground_station.LinearGroundStationMotionSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.ground_station.PiecewiseConstantGroundStationMotionSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.ground_station.CustomGroundStationMotionSettings
   :members:




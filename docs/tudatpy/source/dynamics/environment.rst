.. _environment:

``environment``
===============
Functionalities of environment objects.

This module provides functionalities for environment objects. Specifically, it contains
classes and functions that perform computations related to environment models of natural and artificial bodies.
Much of the functionality in this module concerns classes stored inside :class:`~tudatpy.dynamics.environment.Body` objects, a list of which is in turn
stored in a :class:`~tudatpy.dynamics.environment.SystemOfBodies` object. Note that the classes in this module are rarely created manually, 
but are instead created by the functionality in the :ref:`environment_setup` module. 

Station histories and instrument reference points
------------------------------------------------

Station histories evaluate absolute terrestrial positions in metres and states
in metres and metres per TDB second. Construct the usual IGN trajectory with
``environment_setup.ground_station.ign_station_history`` and attach it with
``station_history_motion``; see :ref:`ground_station` for the complete workflow.
Frame and reference-point construction functions remain in that setup module.

The runtime objects below retain source records, output metadata and explicit
composition order. Returned record lists are independent copies. Source-coordinate
types and the named IGN conversion functions support advanced custom trajectories
and callbacks; the default station-motion settings perform that conversion.

.. currentmodule:: tudatpy.dynamics.environment

.. autoclass:: StationHistoryTrajectory
   :members:

.. autoclass:: TransformedStationHistoryTrajectory
   :members:

.. autoclass:: StationReferencePointTrajectory
   :members:

.. autoclass:: StationReferencePointOffsetHistory
   :members:

.. autoclass:: StationHistoryMetadata
   :members:

.. autoclass:: StationHistoryTime
   :members:

.. autoclass:: SinexEquipmentHeaderTimeContext
   :members:

.. autoclass:: StationReferencePointOffset
   :members:

.. autofunction:: ign_nominal_utc_station_time_from_tdb

.. autofunction:: ign_nominal_utc_station_time_from_tai

.. autofunction:: ign_nominal_utc_station_time

.. autofunction:: sinex_civil_utc_to_tai

.. autofunction:: station_reference_point_identity












Functions
---------
.. currentmodule:: tudatpy.dynamics.environment

.. autosummary::

   transform_to_inertial_orientation
   
   save_vehicle_mesh_to_file

   fit_tle_to_cartesian_state_history



.. autofunction:: tudatpy.dynamics.environment.save_vehicle_mesh_to_file

.. autofunction:: tudatpy.dynamics.environment.transform_to_inertial_orientation

.. autofunction:: tudatpy.dynamics.environment.fit_tle_to_cartesian_state_history


Enums
-----

.. currentmodule:: tudatpy.dynamics.environment

.. autosummary::

   FrequencyGapHandling


.. autoclass:: tudatpy.dynamics.environment.FrequencyGapHandling
   :members:


Classes
-------
.. currentmodule:: tudatpy.dynamics.environment

.. autosummary::


   Ephemeris

   TimeEphemeris

   Tle

   TleFitResult

   RotationalEphemeris

   GcrsToItrsRotationModel

   EarthOrientationAnglesCalculator

   GravityFieldModel

   GravityFieldVariationModel
   
   SphericalHarmonicsGravityField

   TimeDependentSphericalHarmonicsGravityField
   
   BodyShapeModel
   
   RigidBodyProperties
      
   AtmosphereModel
   
   RadiationSourceModel

   AerodynamicCoefficientInterface

   HypersonicLocalInclinationAnalysis

   TransmittingFrequencyCalculator

   ConstantTransmittingFrequencyCalculator

   PiecewiseLinearFrequencyInterpolator

   PointingAnglesCalculator
   
   GroundStation
   
   GroundStationState

   FlightConditions

   AtmosphericFlightConditions

   AerodynamicAngleCalculator
   
   VehicleSystems

   EngineModel

   Body

   SystemOfBodies

   
.. autoclass:: tudatpy.dynamics.environment.Ephemeris
   :members:

.. autoclass:: tudatpy.dynamics.environment.TimeEphemeris
   :members:

.. autoclass:: tudatpy.dynamics.environment.Tle
   :members:
   :special-members: __init__

.. autoclass:: tudatpy.dynamics.environment.TleFitResult
   :members:

.. autoclass:: tudatpy.dynamics.environment.RotationalEphemeris
   :members:

.. autoclass:: tudatpy.dynamics.environment.GcrsToItrsRotationModel
   :members:

.. autoclass:: tudatpy.dynamics.environment.EarthOrientationAnglesCalculator
   :members:

.. autoclass:: tudatpy.dynamics.environment.GravityFieldModel
   :members:

.. autoclass:: tudatpy.dynamics.environment.GravityFieldVariationModel
   :members:

.. autoclass:: tudatpy.dynamics.environment.SphericalHarmonicsGravityField
   :members:
      
.. autoclass:: tudatpy.dynamics.environment.TimeDependentSphericalHarmonicsGravityField
   :members:

.. autoclass:: tudatpy.dynamics.environment.BodyShapeModel
   :members:
  
.. autoclass:: tudatpy.dynamics.environment.RigidBodyProperties
   :members: 

.. autoclass:: tudatpy.dynamics.environment.AtmosphereModel
   :members:     

.. autoclass:: tudatpy.dynamics.environment.RadiationSourceModel
   :members:   

.. autoclass:: tudatpy.dynamics.environment.AerodynamicCoefficientInterface
   :members:

.. autoclass:: tudatpy.dynamics.environment.HypersonicLocalInclinationAnalysis
   :members:
   :special-members: __init__

.. autoclass:: tudatpy.dynamics.environment.TransmittingFrequencyCalculator
   :members:

.. autoclass:: tudatpy.dynamics.environment.ConstantTransmittingFrequencyCalculator
   :members:

.. autoclass:: tudatpy.dynamics.environment.PiecewiseLinearFrequencyInterpolator
   :members:

.. autoclass:: tudatpy.dynamics.environment.PointingAnglesCalculator
   :members:

.. autoclass:: tudatpy.dynamics.environment.GroundStation
   :members:
      
.. autoclass:: tudatpy.dynamics.environment.GroundStationState
   :members:

.. autoclass:: tudatpy.dynamics.environment.FlightConditions
   :members:

.. autoclass:: tudatpy.dynamics.environment.AtmosphericFlightConditions
   :members:

.. autoclass:: tudatpy.dynamics.environment.AerodynamicAngleCalculator
   :members:
   
.. autoclass:: tudatpy.dynamics.environment.VehicleSystems
   :members:
   
.. autoclass:: tudatpy.dynamics.environment.EngineModel
   :members:

.. autoclass:: tudatpy.dynamics.environment.Body
   :members:

.. autoclass:: tudatpy.dynamics.environment.SystemOfBodies
   :members:
   :exclude-members: get_body

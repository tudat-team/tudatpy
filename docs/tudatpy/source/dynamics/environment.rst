.. _environment:

``environment``
===============
Functionalities of environment objects.

This module provides functionalities for environment objects. Specifically, it contains
classes and functions that perform computations related to environment models of natural and artificial bodies.
Much of the functionality in this module concerns classes stored inside :class:`~tudatpy.dynamics.environment.Body` objects, a list of which is in turn
stored in a :class:`~tudatpy.dynamics.environment.SystemOfBodies` object. Note that the classes in this module are rarely created manually, 
but are instead created by the functionality in the :ref:`environment_setup` module. 













Functions
---------
.. currentmodule:: tudatpy.dynamics.environment

.. autosummary::

   transform_to_inertial_orientation
   
   save_vehicle_mesh_to_file

   fit_tle_to_cartesian_state_history

   ign_nominal_utc_station_time_from_tdb

   ign_nominal_utc_station_time_from_tai

   ign_nominal_utc_station_time

   sinex_civil_utc_to_tai

   station_reference_point_identity



.. autofunction:: tudatpy.dynamics.environment.save_vehicle_mesh_to_file

.. autofunction:: tudatpy.dynamics.environment.transform_to_inertial_orientation

.. autofunction:: tudatpy.dynamics.environment.fit_tle_to_cartesian_state_history

.. autofunction:: tudatpy.dynamics.environment.ign_nominal_utc_station_time_from_tdb

.. autofunction:: tudatpy.dynamics.environment.ign_nominal_utc_station_time_from_tai

.. autofunction:: tudatpy.dynamics.environment.ign_nominal_utc_station_time

.. autofunction:: tudatpy.dynamics.environment.sinex_civil_utc_to_tai

.. autofunction:: tudatpy.dynamics.environment.station_reference_point_identity


Enums
-----

.. currentmodule:: tudatpy.dynamics.environment

.. autosummary::

   FrequencyGapHandling

   StationReferencePointKind

   StationPostSeismicFunction

   StationLocalVertical

   StationReferencePointOffsetBasis


.. autoclass:: tudatpy.dynamics.environment.FrequencyGapHandling
   :members:

.. autoclass:: tudatpy.dynamics.environment.StationReferencePointKind
   :members:

.. autoclass:: tudatpy.dynamics.environment.StationPostSeismicFunction
   :members:

.. autoclass:: tudatpy.dynamics.environment.StationLocalVertical
   :members:

.. autoclass:: tudatpy.dynamics.environment.StationReferencePointOffsetBasis
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

   StationCoordinateHistory

   StationReferencePointOffsetHistory

   StationReferencePointOffset

   StationReferencePointIdentity

   StationReferencePointGeometry

   StationHistoryMetadata

   StationHistoryTime

   StationHistoryInterval

   StationHistoryLinearSolution

   StationPostSeismicTerm

   StationSeasonalTerm

   SinexEquipmentHeaderTimeContext

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

.. autoclass:: tudatpy.dynamics.environment.StationCoordinateHistory
   :members:

.. autoclass:: tudatpy.dynamics.environment.StationReferencePointOffsetHistory
   :members:

.. autoclass:: tudatpy.dynamics.environment.StationReferencePointOffset
   :members:

.. autoclass:: tudatpy.dynamics.environment.StationReferencePointIdentity
   :members:

.. autoclass:: tudatpy.dynamics.environment.StationReferencePointGeometry
   :members:

.. autoclass:: tudatpy.dynamics.environment.StationHistoryMetadata
   :members:

.. autoclass:: tudatpy.dynamics.environment.StationHistoryTime
   :members:

.. autoclass:: tudatpy.dynamics.environment.StationHistoryInterval
   :members:

.. autoclass:: tudatpy.dynamics.environment.StationHistoryLinearSolution
   :members:

.. autoclass:: tudatpy.dynamics.environment.StationPostSeismicTerm
   :members:

.. autoclass:: tudatpy.dynamics.environment.StationSeasonalTerm
   :members:

.. autoclass:: tudatpy.dynamics.environment.SinexEquipmentHeaderTimeContext
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

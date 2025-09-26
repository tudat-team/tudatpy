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



.. autofunction:: tudatpy.dynamics.environment.save_vehicle_mesh_to_file

.. autofunction:: tudatpy.dynamics.environment.transform_to_inertial_orientation


Classes
-------
.. currentmodule:: tudatpy.dynamics.environment

.. autosummary::


   Ephemeris
   
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




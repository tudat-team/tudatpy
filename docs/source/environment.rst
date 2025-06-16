.. _environment:

``environment``
===============
Functionalities of environment objects.

This module provides functionalities for environment objects. Specifically, it contains
classes and functions that perform computations related to environment models of natural and artificial bodies.
Much of the functionality in this module concerns classes stored inside :class:`~tudatpy.numerical_simulation.environment.Body` objects, a list of which is in turn
stored in a :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies` object. Note that the classes in this module are rarely created manually, 
but are instead created by the functionality in the :ref:`environment_setup` module. 












Functions
---------
.. currentmodule:: tudatpy.numerical_simulation.environment

.. autosummary::

   transform_to_inertial_orientation
   
   save_vehicle_mesh_to_file



.. autofunction:: tudatpy.numerical_simulation.environment.save_vehicle_mesh_to_file

.. autofunction:: tudatpy.numerical_simulation.environment.transform_to_inertial_orientation



Enumerations
------------
.. currentmodule:: tudatpy.numerical_simulation.environment

.. autosummary::

   AerodynamicsReferenceFrames

   AerodynamicCoefficientFrames

   AerodynamicCoefficientsIndependentVariables
   
   AtmosphericCompositionSpecies


.. autoclass:: tudatpy.numerical_simulation.environment.AerodynamicsReferenceFrameAngles
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.AerodynamicsReferenceFrames
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.AerodynamicCoefficientFrames
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.AerodynamicCoefficientsIndependentVariables
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.AtmosphericCompositionSpecies
   :members:


Classes
-------
.. currentmodule:: tudatpy.numerical_simulation.environment

.. autosummary::


   Ephemeris
   
   RotationalEphemeris

   GcrsToItrsRotationModel

   EarthOrientationAnglesCalculator

   GravityFieldModel
   
   SphericalHarmonicsGravityField
   
   BodyShapeModel
   
   RigidBodyProperties
      
   AtmosphereModel
   
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



   
.. autoclass:: tudatpy.numerical_simulation.environment.Ephemeris
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.RotationalEphemeris
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.GcrsToItrsRotationModel
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.EarthOrientationAnglesCalculator
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.GravityFieldModel
   :members:
      
.. autoclass:: tudatpy.numerical_simulation.environment.SphericalHarmonicsGravityField
   :members:
      
.. autoclass:: tudatpy.numerical_simulation.environment.BodyShapeModel
   :members:
  
.. autoclass:: tudatpy.numerical_simulation.environment.RigidBodyProperties
   :members: 

.. autoclass:: tudatpy.numerical_simulation.environment.AtmosphereModel
   :members:     
   
.. autoclass:: tudatpy.numerical_simulation.environment.AerodynamicCoefficientInterface
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.HypersonicLocalInclinationAnalysis
   :members:
   :special-members: __init__

.. autoclass:: tudatpy.numerical_simulation.environment.GroundStation
   :members:
      
.. autoclass:: tudatpy.numerical_simulation.environment.GroundStationState
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.FlightConditions
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.AtmosphericFlightConditions
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.AerodynamicAngleCalculator
   :members:
   
.. autoclass:: tudatpy.numerical_simulation.environment.VehicleSystems
   :members:
   
.. autoclass:: tudatpy.numerical_simulation.environment.EngineModel
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.Body
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.SystemOfBodies
   :members:




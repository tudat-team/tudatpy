``environment``
===============
Functionalities of environment objects.

This module provides functionalities for environment objects. Specifically, it contains
classes and functions that perform computations related to environment models of natural and artificial bodies.
Much of the functionality in this module concerns classes stored inside :class:`~tudatpy.numerical_simulation.environment.Body` objects, a list of which is in turn
stored in a :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies` object. Note that the classes in this module are rarely created manually, 
but are instead created by the functionality in the :ref:`\`\`environment_setup\`\`` module. 












Functions
---------
.. currentmodule:: tudatpy.numerical_simulation.environment

.. autosummary::

   save_vehicle_mesh_to_file



.. autofunction:: tudatpy.numerical_simulation.environment.save_vehicle_mesh_to_file




Enumerations
------------
.. currentmodule:: tudatpy.numerical_simulation.environment

.. autosummary::

   AerodynamicsReferenceFrames

   AerodynamicCoefficientFrames

   AerodynamicCoefficientsIndependentVariables



.. autoclass:: tudatpy.numerical_simulation.environment.AerodynamicsReferenceFrames
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.AerodynamicCoefficientFrames
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.AerodynamicCoefficientsIndependentVariables
   :members:




Classes
-------
.. currentmodule:: tudatpy.numerical_simulation.environment

.. autosummary::

   AerodynamicCoefficientInterface

   HypersonicLocalInclinationAnalysis

   FlightConditions

   AtmosphericFlightConditions

   AerodynamicAngleCalculator

   RotationalEphemeris

   VehicleSystems

   Ephemeris

   EngineModel

   Body

   SystemOfBodies



.. autoclass:: tudatpy.numerical_simulation.environment.AerodynamicCoefficientInterface
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.HypersonicLocalInclinationAnalysis
   :members:
   :special-members: __init__



.. autoclass:: tudatpy.numerical_simulation.environment.FlightConditions
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.AtmosphericFlightConditions
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.AerodynamicAngleCalculator
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.RotationalEphemeris
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.VehicleSystems
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.Ephemeris
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.EngineModel
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.Body
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment.SystemOfBodies
   :members:




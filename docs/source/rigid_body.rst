``rigid_body``
==============
This module contains a set of factory functions for setting up the
models for a rigid body in Tudat. Specifically, this module defines the
mass, center of mass and inertia tensor of a body. These models **do not**
define a gravity field of a body. However, when a body is endowed with a 
gravity field, a compatible rigid body is automatically created. 












Functions
---------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.rigid_body

.. autosummary::

   constant_rigid_body_properties

   custom_time_dependent_rigid_body_properties

   custom_mass_dependent_rigid_body_properties



.. autofunction:: tudatpy.numerical_simulation.environment_setup.rigid_body.constant_rigid_body_properties

.. autofunction:: tudatpy.numerical_simulation.environment_setup.rigid_body.custom_time_dependent_rigid_body_properties

.. autofunction:: tudatpy.numerical_simulation.environment_setup.rigid_body.custom_mass_dependent_rigid_body_properties






Classes
-------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.rigid_body

.. autosummary::

   RigidBodyPropertiesSettings



.. autoclass:: tudatpy.numerical_simulation.environment_setup.rigid_body.RigidBodyPropertiesSettings
   :members:




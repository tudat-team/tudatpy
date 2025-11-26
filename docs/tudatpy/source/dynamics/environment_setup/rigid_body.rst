.. _rigid_body:

``rigid_body``
==============
This module contains a set of factory functions for setting up the
models for so-called 'rigid body properties' in Tudat. Specifically, rigid body properties define the
mass, center of mass and inertia tensor of a body.

The main interfaces with Tudat is the :attr:`~tudatpy.dynamics.environment_setup.BodySettings.rigid_body_settings`
attribute  (of type :class:`~tudatpy.dynamics.environment_setup.rigid_body.RigidBodyPropertiesSettings`) of the body settings, which defines settings for the rigid-body properties of a body.
**The functions in this submodule are used to create these settings objects.** When creating a body (typically using the
:func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies` function), an object of type
:class:`~tudatpy.dynamics.environment.RigidBodyProperties` (or a derived class) is created
and added to the associated :class:`~tudatpy.dynamics.environment.Body` object based on the settings object, which can
be retrieved using the :attr:`~tudatpy.dynamics.environment.Body.rigid_body_properties` attribute.

The models defined here **do not** define a gravity field of a body. However, when a body is endowed with a
gravity field (see :ref:`gravity_field`), compatible rigid body settings are automatically created, mapping the gravitational
paramater of the gravity field to a mass here. For a :func:`~tudatpy.dynamics.environment_setup.gravity_field.spherical_harmonic`
gravity field, the degree-one coefficients are used to set the body center of mass.
In addition, when creating a spherica harmonic gravity field, and specifying a :attr:`~tudatpy.dynamics.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings.scaled_mean_moment_of_inertia`, an inertia tensor is created and defined in the rigid-body settings.

Functions
---------
.. currentmodule:: tudatpy.dynamics.environment_setup.rigid_body

.. autosummary::

   constant_rigid_body_properties

   custom_time_dependent_rigid_body_properties

   custom_mass_dependent_rigid_body_properties



.. autofunction:: tudatpy.dynamics.environment_setup.rigid_body.constant_rigid_body_properties

.. autofunction:: tudatpy.dynamics.environment_setup.rigid_body.custom_time_dependent_rigid_body_properties

.. autofunction:: tudatpy.dynamics.environment_setup.rigid_body.custom_mass_dependent_rigid_body_properties






Classes
-------
.. currentmodule:: tudatpy.dynamics.environment_setup.rigid_body

.. autosummary::

   RigidBodyPropertiesSettings



.. autoclass:: tudatpy.dynamics.environment_setup.rigid_body.RigidBodyPropertiesSettings
   :members:




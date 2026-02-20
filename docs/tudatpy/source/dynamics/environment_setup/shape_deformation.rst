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

The shape deformation models are used for global shape variations of a celesial body, which are used (for instance)
for the high-accuracy modelling of the positions of ground stations on a body. Station-specific models for position variations
of a station (such as plate motion) are set in the ground stations themselves (see :ref:`ground_station`).


Functions
---------
.. currentmodule:: tudatpy.dynamics.environment_setup.shape_deformation

.. autosummary::

   basic_solid_body_tidal

   degree_two_basic_solid_body_tidal

   iers_2010_solid_body_tidal

   ocean_tide

   pole_tide


.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.basic_solid_body_tidal

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.degree_two_basic_solid_body_tidal

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.iers_2010_solid_body_tidal

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.ocean_tide

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.pole_tide





Classes
-------
.. currentmodule:: tudatpy.dynamics.environment_setup.shape_deformation

.. autosummary::

   BodyDeformationSettings

   BasicSolidBodyDeformationSettings



.. autoclass:: tudatpy.dynamics.environment_setup.shape_deformation.BodyDeformationSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.shape_deformation.BasicSolidBodyDeformationSettings
   :members:




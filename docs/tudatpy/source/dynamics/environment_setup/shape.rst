.. _shape:

``shape``
=========

This module contains a set of factory functions for setting up the
physical shape of natural bodies in an environment.

The main interfaces with Tudat is the :attr:`~tudatpy.dynamics.environment_setup.BodySettings.shape_settings`
attribute  (of type :class:`~tudatpy.dynamics.environment_setup.shape.BodyShapeSettings`) of the body settings, which defines settings for the shape model of a body.
**The functions in this submodule are used to create these settings objects.** When creating a body (typically using the
:func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies` function), an object of type
:class:`~tudatpy.dynamics.environment.BodyShapeModel` (or a derived class) is created
and added to the associated :class:`~tudatpy.dynamics.environment.Body` object based on the settings object, which can
be retrieved using the :attr:`~tudatpy.dynamics.environment.Body.shape_model` attribute.

The following code block gives an overview of the steps to define, create, and extract a shape model, for the specific example of an oblate spheroid shape model
  with :math:`R_{e}=3396.2` km equatorial radius, and flattening :math:`f=0.00589`

.. code-block:: python

  from tudatpy.dynamics import environment_setup
  from tudatpy.astro import time_representation

  # Create body settings
  body_settings =  environment_setup.get_default_body_settings( ... ) # Typical way to instantiate body settings

  # Modify shape model settings (base class type BodyShapeSettings)
  body_settings.get( 'Mars' ).shape_settings = environment_setup.shape.oblate_spherical(
      equatorial_radius = 3396.2E3,
      flattening = 0.00589 )

  # Create bodies
  bodies = environment_setup.create_system_of_bodies(body_settings)

  # Extract shape model (base class type ShapeModel) from Earth
  mars_shape_model = bodies.get( 'Mars' ).shape_model


Note that, in Tudat, the physical shape of a body and its gravity field are completely independent. In addition,
the functionality in this module is for defining shapes of natural bodies. Shapes of spacecraft (e.g. macromodels)
are defined through functionality in the :func:`vehicle_systems` module.





Functions
---------
.. currentmodule:: tudatpy.dynamics.environment_setup.shape

.. autosummary::

   spherical

   spherical_spice

   oblate_spherical

   polyhedron

   hybrid



.. autofunction:: tudatpy.dynamics.environment_setup.shape.spherical

.. autofunction:: tudatpy.dynamics.environment_setup.shape.spherical_spice

.. autofunction:: tudatpy.dynamics.environment_setup.shape.oblate_spherical

.. autofunction:: tudatpy.dynamics.environment_setup.shape.polyhedron

.. autofunction:: tudatpy.dynamics.environment_setup.shape.hybrid






Classes
-------
.. currentmodule:: tudatpy.dynamics.environment_setup.shape

.. autosummary::

   BodyShapeSettings

   SphericalBodyShapeSettings

   OblateSphericalBodyShapeSettings

   PolyhedronBodyShapeSettings

   HybridBodyShapeSettings



.. autoclass:: tudatpy.dynamics.environment_setup.shape.BodyShapeSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.shape.SphericalBodyShapeSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.shape.OblateSphericalBodyShapeSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.shape.PolyhedronBodyShapeSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.shape.HybridBodyShapeSettings
   :members:




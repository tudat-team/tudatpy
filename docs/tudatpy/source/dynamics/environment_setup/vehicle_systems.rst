.. _vehicle_systems:

``vehicle_systems``
===================
This module contains a set of factory functions for setting up physical and system properties of a vehicle.
For various high-accuracy models of non-conservative spacecraft dynamics, a so-called macromodel is required which defines
the external shape of the vehicle. This macromodel is typically defined by a set of panels, with each panel assigned
specific properties of how it interacts with the environment. At present, the spacecraft macromodel in Tudat is
used for the calculation of a panelled radiation pressure acceleration and the
calculation of aerodynamic coefficients in both rarefied and hypersonic flow using gas-surface interaction models (GSMIs).

The main interfaces with Tudat for a vehicle shape model is the :attr:`~tudatpy.dynamics.environment_setup.BodySettings.vehicle_shape_settings`
attribute  (of type :class:`~tudatpy.dynamics.environment_setup.vehicle_systems.FullPanelledBodySettings`) of the body settings, which defines settings for the macromodel of a body.
**The functions in this submodule are used to create these settings objects.** When creating a body (typically using the
:func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies` function), an object of type
:class:`~tudatpy.dynamics.environment.VehicleSystems` is created
and added to the associated :class:`~tudatpy.dynamics.environment.Body` object based on the settings object, which can
be retrieved using the :attr:`~tudatpy.dynamics.environment.Body.system_models` attribute. The ``VehicleSystems`` contains (among various other models) the macromodel properties of the body.

The current panels in Tudat allow a list of panels to be defined, with the minimal geometrical properties of panel :math:`i` defined by the
surface normal vector :math:`\hat{\mathbf{n}}_{i}` and the surface area :math:`A_{i}`. Another way to completely define the geometry of a panel in space is by defining its three vertices 
(triangles allow for speed-ups at the algorithm level), a rotation frame and its origin in the body-fixed frame, allowing for the computation of self-shadowing (SSH) and panelled aerodynamic coefficients. 
Due to the complexity of manually defining each triangular panel, this option is available only by loading a full macromodel.

The panel surface normal may be defined in either the body-fixed frame :math:`\mathcal{B}` of the vehicle, or to a 'vehicle-part-fixed frame'
:math:`\mathcal{F}_{j}`. A 'vehicle part' is defined as a part of the vehicle that can move/rotate w.r.t. the body-fixed frame of the
spacecraft. Typical examples are the solar arrays and an movable antenna.

The panel surface normal (in either the body frame or the part frame), may be defined by the
:func:`~tudatpy.dynamics.environment_setup.vehicle_systems.frame_fixed_panel_geometry`,
:func:`~tudatpy.dynamics.environment_setup.vehicle_systems.time_varying_panel_geometry` or
:func:`~tudatpy.dynamics.environment_setup.vehicle_systems.body_tracking_panel_geometry` functions,
where the latter is used to ensure that a panel normal automatically points to/away from another bodY (e.g. the Sun for solar panels).

A full panel is created by defining its geometry, and models for its interaction with the environment (currently limited to
a reflection law to compute the influence of radiation pressure) using the
:func:`~tudatpy.dynamics.environment_setup.vehicle_systems.body_panel_settings` function.

The vehicle macromodel, and the rotation models from the body-fixed frame to the (optional) part-fixed frames are defined by
using the :func:`~tudatpy.dynamics.environment_setup.vehicle_systems.full_panelled_body_settings` function, and
assigned to the :attr:`~tudatpy.dynamics.environment_setup.BodySettings.vehicle_shape_settings` attribute of the :class:`~tudatpy.dynamics.environment_setup.BodySettings` class. Large macromodels can be loaded from CAD files in DAE format using :func:`~tudatpy.dynamics.environment_setup.vehicle_systems.body_panel_settings_list_from_dae`.
When a full macromodel is not available to the user, a 'box-wing' model may also be used, which creates the macromodel
based on user settings, using the :func:`~tudatpy.dynamics.environment_setup.vehicle_systems.box_wing_panelled_body_settings` function.












Functions
---------
.. currentmodule:: tudatpy.dynamics.environment_setup.vehicle_systems

.. autosummary::

   frame_fixed_panel_geometry

   time_varying_panel_geometry

   body_tracking_panel_geometry

   material_properties

   body_panel_settings

   full_panelled_body_settings

   box_wing_panelled_body_settings

   body_panel_settings_list_from_dae



.. autofunction:: tudatpy.dynamics.environment_setup.vehicle_systems.frame_fixed_panel_geometry

.. autofunction:: tudatpy.dynamics.environment_setup.vehicle_systems.time_varying_panel_geometry

.. autofunction:: tudatpy.dynamics.environment_setup.vehicle_systems.body_tracking_panel_geometry

.. autofunction:: tudatpy.dynamics.environment_setup.vehicle_systems.material_properties

.. autofunction:: tudatpy.dynamics.environment_setup.vehicle_systems.body_panel_settings

.. autofunction:: tudatpy.dynamics.environment_setup.vehicle_systems.full_panelled_body_settings

.. autofunction:: tudatpy.dynamics.environment_setup.vehicle_systems.box_wing_panelled_body_settings

.. autofunction:: tudatpy.dynamics.environment_setup.vehicle_systems.body_panel_settings_list_from_dae






Classes
-------
.. currentmodule:: tudatpy.dynamics.environment_setup.vehicle_systems

.. autosummary::

   MaterialProperties

   BodyPanelGeometrySettings

   FrameFixedBodyPanelGeometrySettings

   FrameVariableBodyPanelGeometrySettings

   BodyPanelSettings

   FullPanelledBodySettings


.. autoclass:: tudatpy.dynamics.environment_setup.vehicle_systems.MaterialProperties
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.vehicle_systems.BodyPanelGeometrySettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.vehicle_systems.FrameFixedBodyPanelGeometrySettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.vehicle_systems.FrameVariableBodyPanelGeometrySettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.vehicle_systems.BodyPanelSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.vehicle_systems.FullPanelledBodySettings
   :members:




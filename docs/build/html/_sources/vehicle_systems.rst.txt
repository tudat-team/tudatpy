``vehicle_systems``
===================
This module contains a set of factory functions for setting up physical and system properties of a vehicle.
For various high-accuracy models of non-conservative spacecraft dynamics, a so-called macromodel is required which defines
the external shape of the vehicle. This macromodel is typically defined by a set of panels, with each panel assigned
specific properties of how it interacts with the environment. At present, the spacecraft macromodel in Tudat is only
used for the calculation of a panelled radiation pressure acceleration, but future updates will also use it for the
calculation of aerodynamic coefficients in both rarefied and hypersonic flow.

The current panels in Tudat allow a list of panels to be defined, with the geometrical properties of panel :math:`i` defined by the
surface normal vector :math:`\hat{\mathbf{n}}_{i}` and the surface area :math:`A_{i}`. Note that, since the panel shape or
location is not yet defined, computing torques due to surface forces, or incorporating shadowing into the panel
force calculatuion, is not yet supported.

The panel surface normal may be defined in either the body-fixed frame :math:`\mathcal{B}` of the vehicle, or to a 'vehicle-part-fixed frame'
:math:`\mathcal{F}_{j}`. A 'vehicle part' is defined as a part of the vehicle that can move/rotate w.r.t. the body-fixed frame of the
spacecraft. Typical examples are the solar arrays and an movable antenna.

The panel surface normal (in either the body frame or the part frame), may be defined by the
:func:`~tudatpy.numerical_simulation.environment_setup.vehicle_systems.frame_fixed_panel_geometry`,
:func:`~tudatpy.numerical_simulation.environment_setup.vehicle_systems.time_varying_panel_geometry` or
:func:`~tudatpy.numerical_simulation.environment_setup.vehicle_systems.body_tracking_panel_geometry` functions,
where the latter is used to ensure that a panel normal automatically points to/away from another bodY (e.g. the Sun for solar panels).

A full panel is created by defining its geometry, and models for its interaction with the environment (currently limited to
a reflection law to compute the influence of radiation pressure) using the
:func:`~tudatpy.numerical_simulation.environment_setup.vehicle_systems.body_panel_settings` function.

The vehicle macromodel, and the rotation models from the body-fixed frame to the (optional) part-fixed frames are defined by
using the :func:`~tudatpy.numerical_simulation.environment_setup.vehicle_systems.full_panelled_body_settings` function, and
assigned to the ``vehicle_shape_settings`` attribute of the :class:`~tudatpy.numerical_simulation.environment_setup.BodySettings` class.
When a full macromodel is not available to the user, a 'box-wing' model may also be used, which creates the macromodel
bassed on user settings, using the :func:`~tudatpy.numerical_simulation.environment_setup.vehicle_systems.box_wing_panelled_body_settings` function.












Functions
---------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.vehicle_systems

.. autosummary::

   frame_fixed_panel_geometry

   time_varying_panel_geometry

   body_tracking_panel_geometry

   body_panel_settings

   full_panelled_body_settings

   box_wing_panelled_body_settings



.. autofunction:: tudatpy.numerical_simulation.environment_setup.vehicle_systems.frame_fixed_panel_geometry

.. autofunction:: tudatpy.numerical_simulation.environment_setup.vehicle_systems.time_varying_panel_geometry

.. autofunction:: tudatpy.numerical_simulation.environment_setup.vehicle_systems.body_tracking_panel_geometry

.. autofunction:: tudatpy.numerical_simulation.environment_setup.vehicle_systems.body_panel_settings

.. autofunction:: tudatpy.numerical_simulation.environment_setup.vehicle_systems.full_panelled_body_settings

.. autofunction:: tudatpy.numerical_simulation.environment_setup.vehicle_systems.box_wing_panelled_body_settings






Classes
-------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.vehicle_systems

.. autosummary::

   BodyPanelGeometrySettings

   BodyPanelSettings

   FullPanelledBodySettings



.. autoclass:: tudatpy.numerical_simulation.environment_setup.vehicle_systems.BodyPanelGeometrySettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.vehicle_systems.BodyPanelSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.vehicle_systems.FullPanelledBodySettings
   :members:




.. _atmosphere:

``atmosphere``
==============
This module contains a set of factory functions for setting up the
atmosphere models of celestial bodies in an environment.

The main interfaces with Tudat is the :attr:`~tudatpy.dynamics.environment_setup.BodySettings.atmosphere_settings`
attribute (of type :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings`) of the body settings, which defines settings for the atmosphere of a body.
**The functions in this submodule are used to create such settings objects**. When creating a body (typically using the
:func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies` function), an object of type
:class:`~tudatpy.dynamics.environment.AtmosphereModel` (or a derived class) is created
and added to the associated :class:`~tudatpy.dynamics.environment.Body` object based on the settings object, which can
be retrieved using the :attr:`~tudatpy.dynamics.environment.Body.atmosphere_model` attribute.

As a property of the atmosphere model, Tudat allows a wind model to be defined. In the absence of a wind model, the
atmosphere is assumed to co-rotate with the central body. Specifically, the wind velocity in the central-body fixed frame
is exactly :math:`\mathbf{0}` by default. Settings for a wind model to deviate from this are added to the atmosphere settings
through the :attr:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings.wind_settings` attribute (of type :class:`~tudatpy.dynamics.environment_setup.atmosphere.WindModelSettings`).
Several functions in this submodule allow an object of this type to be created.

The atmospheric properties of a body are used in a number of Tudat models, but the primary impact on a numerical
propagation/estimation is through the impact of atmospheric density :math:`\rho`, which is used by the
:func:`~tudatpy.dynamics.propagation_setup.acceleration.aerodynamic` acceleration model and
:func:`~tudatpy.dynamics.propagation_setup.torque.aerodynamic` torque model.







Functions
---------
.. currentmodule:: tudatpy.dynamics.environment_setup.atmosphere

.. autosummary::

   constant_wind_model
   custom_wind_model
   exponential_predefined
   exponential
   nrlmsise00
   tabulated
   us76
   mars_dtm
   custom_constant_temperature
   custom_four_dimensional_constant_temperature
   scaled_by_constant
   scaled_by_function


.. autofunction:: tudatpy.dynamics.environment_setup.atmosphere.constant_wind_model

.. autofunction:: tudatpy.dynamics.environment_setup.atmosphere.custom_wind_model

.. autofunction:: tudatpy.dynamics.environment_setup.atmosphere.exponential_predefined

.. autofunction:: tudatpy.dynamics.environment_setup.atmosphere.exponential

.. autofunction:: tudatpy.dynamics.environment_setup.atmosphere.nrlmsise00

.. autofunction:: tudatpy.dynamics.environment_setup.atmosphere.tabulated

.. autofunction:: tudatpy.dynamics.environment_setup.atmosphere.us76

.. autofunction:: tudatpy.dynamics.environment_setup.atmosphere.mars_dtm

.. autofunction:: tudatpy.dynamics.environment_setup.atmosphere.custom_constant_temperature

.. autofunction:: tudatpy.dynamics.environment_setup.atmosphere.custom_four_dimensional_constant_temperature

.. autofunction:: tudatpy.dynamics.environment_setup.atmosphere.scaled_by_constant

.. autofunction:: tudatpy.dynamics.environment_setup.atmosphere.scaled_by_function


Classes
-------
.. currentmodule:: tudatpy.dynamics.environment_setup.atmosphere

.. autosummary::

   WindModelSettings
   ConstantWindModelSettings
   CustomWindModelSettings
   AtmosphereSettings
   ExponentialAtmosphereSettings
   CustomConstantTemperatureAtmosphereSettings
   ScaledAtmosphereSettings
   NRLMSISE00Input
   NRLMSISE00Atmosphere


Wind Model Settings
~~~~~~~~~~~~~~~~~~~

.. autoclass:: tudatpy.dynamics.environment_setup.atmosphere.WindModelSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.atmosphere.ConstantWindModelSettings
   :members:
   :show-inheritance:

.. autoclass:: tudatpy.dynamics.environment_setup.atmosphere.CustomWindModelSettings
   :members:
   :show-inheritance:


Atmosphere Settings
~~~~~~~~~~~~~~~~~~~

.. autoclass:: tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.atmosphere.ExponentialAtmosphereSettings
   :members:
   :show-inheritance:

.. autoclass:: tudatpy.dynamics.environment_setup.atmosphere.CustomConstantTemperatureAtmosphereSettings
   :members:
   :show-inheritance:

.. autoclass:: tudatpy.dynamics.environment_setup.atmosphere.ScaledAtmosphereSettings
   :members:
   :show-inheritance:


Atmosphere Model Classes
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: tudatpy.dynamics.environment_setup.atmosphere.NRLMSISE00Input
   :members:
   :special-members: __init__

.. autoclass:: tudatpy.dynamics.environment_setup.atmosphere.NRLMSISE00Atmosphere
   :members:
   :special-members: __init__


Enumerations
------------
.. currentmodule:: tudatpy.dynamics.environment_setup.atmosphere

.. autosummary::

   AtmosphereDependentVariables


.. autoclass:: tudatpy.dynamics.environment_setup.atmosphere.AtmosphereDependentVariables
   :members:
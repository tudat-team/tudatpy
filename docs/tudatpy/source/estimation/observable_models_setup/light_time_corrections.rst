.. _light_time_corrections:

``light_time_corrections``
===========================

This module contains a set of factory functions for setting up the
light-time corrections for observation models.

Most functions in this module create objects of type :class:`~tudatpy.estimation.biases.LightTimeCorrectionSettings`,
which define settings for a type of light-time correction. The main interface with Tudat is that these objects are used
as input to the observation model functions created in the :ref:`model_settings` module.

When not applying any light-time corrections, a signal is assumed to travel in a straight line (in Euclidean space) with
the speed of light :math:`c`. The ligh-time corrections defined through this module are used to compute corrections
:math:`\Delta t` to the light time, such that:

.. math::

   \frac{||\mathbf{r}_{1}(t_{1}) - \mathbf{r}_{0}(t_{0})||}{c}=\left(t_{1}-t_{0}\right)+\Delta t(t_{0},t_{1};\mathbf{r}_{1}(t_{1}),\mathbf{r}_{0}(t_{0}))

for a signal transmitted from link end :math:`0` at :math:`t_{0}` and received by link end :math:`1` at :math:`t_{1}`. More details on the
link with the observation model, and the manner in which this equation is solved, is given on the `user guide <https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-model-setup.html>`_ .

Functions
---------
.. currentmodule:: tudatpy.estimation.observable_models_setup.light_time_corrections

.. autosummary::

   light_time_convergence_settings

   first_order_relativistic_light_time_correction

   approximated_second_order_relativistic_light_time_correction

   dsn_tabulated_tropospheric_light_time_correction

   saastamoinen_tropospheric_light_time_correction

   dsn_tabulated_ionospheric_light_time_correction

   jakowski_ionospheric_light_time_correction

   ionex_ionospheric_light_time_correction

   vmf3_tropospheric_light_time_correction

   inverse_power_series_solar_corona_light_time_correction

   set_vmf_troposphere_data

   set_ionosphere_model_from_ionex

.. autofunction:: tudatpy.estimation.observable_models_setup.light_time_corrections.light_time_convergence_settings

.. autofunction:: tudatpy.estimation.observable_models_setup.light_time_corrections.first_order_relativistic_light_time_correction

.. autofunction:: tudatpy.estimation.observable_models_setup.light_time_corrections.approximated_second_order_relativistic_light_time_correction

.. autofunction:: tudatpy.estimation.observable_models_setup.light_time_corrections.dsn_tabulated_tropospheric_light_time_correction

.. autofunction:: tudatpy.estimation.observable_models_setup.light_time_corrections.saastamoinen_tropospheric_light_time_correction

.. autofunction:: tudatpy.estimation.observable_models_setup.light_time_corrections.dsn_tabulated_ionospheric_light_time_correction

.. autofunction:: tudatpy.estimation.observable_models_setup.light_time_corrections.jakowski_ionospheric_light_time_correction

.. autofunction:: tudatpy.estimation.observable_models_setup.light_time_corrections.ionex_ionospheric_light_time_correction

.. autofunction:: tudatpy.estimation.observable_models_setup.light_time_corrections.vmf3_tropospheric_light_time_correction

.. autofunction:: tudatpy.estimation.observable_models_setup.light_time_corrections.inverse_power_series_solar_corona_light_time_correction

.. autofunction:: tudatpy.estimation.observable_models_setup.light_time_corrections.set_vmf_troposphere_data

.. autofunction:: tudatpy.estimation.observable_models_setup.light_time_corrections.set_ionosphere_model_from_ionex

Enumerations
------------
.. currentmodule:: tudatpy.estimation.observable_models_setup.light_time_corrections

.. autosummary::

   LightTimeFailureHandling

   TroposphericMappingModel

   WaterVaporPartialPressureModel


.. autoclass:: tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeFailureHandling
   :members:

.. autoclass:: tudatpy.estimation.observable_models_setup.light_time_corrections.TroposphericMappingModel
   :members:

.. autoclass:: tudatpy.estimation.observable_models_setup.light_time_corrections.WaterVaporPartialPressureModel
   :members:


Classes
-------
.. currentmodule:: tudatpy.estimation.observable_models_setup.light_time_corrections

.. autosummary::

   LightTimeConvergenceCriteria

   LightTimeCorrectionSettings

   JakowskiVtecCalculator

   GlobalIonosphereModelVtecCalculator

.. autoclass:: tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeConvergenceCriteria
   :members:

.. autoclass:: tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings
   :members:

.. autoclass:: tudatpy.estimation.observable_models_setup.light_time_corrections.JakowskiVtecCalculator
   :members:
   :special-members: __init__

.. autoclass:: tudatpy.estimation.observable_models_setup.light_time_corrections.GlobalIonosphereModelVtecCalculator
   :members:
   :special-members: __init__

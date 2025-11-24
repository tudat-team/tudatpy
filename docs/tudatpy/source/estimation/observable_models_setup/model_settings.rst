.. _model_settings:

``model_settings``
==================

This module contains a set of factory functions for creating settings for observation models (e.g. one-way range,
two-way Doppler, angular position, etc.). The functions here all create an object of type :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings`, from which observation models are created.

Many observation models are based on the transfer of electromagnetic signals (e.g. laser ranging, radio Doppler, etc.),
which require the solution of the light-time solution. Deviations from ideal Euclidean propagation at the speed of light
can be modelled using light-time corrections, for which settings (as objects of type :class:`~tudatpy.estimation.light_time_correction..LightTimeCorrectionSettings`) are created in the :ref:`light_time_corrections` module and provided as input to functions here. Similarly, settings
for observation biases (as objects of type :class:`~tudatpy.estimation.observable_models_setup.biases.ObservationBiasSettings`) are created
in the :ref:`biases` module. Definition of the observation link are created through the :ref:`links` module.

Details on the procedure to create observation models, and the various top-level options for their creatoon, is given on the `user guide <https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-model-setup.html#defining-observation-settings>`_.





Functions
---------
.. currentmodule:: tudatpy.estimation.observable_models_setup.model_settings
.. autosummary::

   one_way_range

   two_way_range

   two_way_range_from_one_way_links

   n_way_range

   n_way_range_from_one_way_links

   dsn_n_way_range

   one_way_doppler_averaged

   two_way_doppler_averaged

   two_way_doppler_averaged_from_one_way_links

   n_way_doppler_averaged

   n_way_doppler_averaged_from_one_way_links

   dsn_n_way_doppler_averaged

   one_way_doppler_instantaneous

   two_way_doppler_instantaneous

   two_way_doppler_instantaneous_from_one_way_links

   angular_position

   relative_angular_position

   cartesian_position

   relative_cartesian_position

   cartesian_velocity

   euler_angles_313


.. autofunction:: tudatpy.estimation.observable_models_setup.model_settings.one_way_range

.. autofunction:: tudatpy.estimation.observable_models_setup.model_settings.two_way_range

.. autofunction:: tudatpy.estimation.observable_models_setup.model_settings.two_way_range_from_one_way_links

.. autofunction:: tudatpy.estimation.observable_models_setup.model_settings.n_way_range

.. autofunction:: tudatpy.estimation.observable_models_setup.model_settings.n_way_range_from_one_way_links

.. autofunction:: tudatpy.estimation.observable_models_setup.model_settings.dsn_n_way_range

.. autofunction:: tudatpy.estimation.observable_models_setup.model_settings.one_way_doppler_averaged

.. autofunction:: tudatpy.estimation.observable_models_setup.model_settings.two_way_doppler_averaged

.. autofunction:: tudatpy.estimation.observable_models_setup.model_settings.two_way_doppler_averaged_from_one_way_links

.. autofunction:: tudatpy.estimation.observable_models_setup.model_settings.n_way_doppler_averaged

.. autofunction:: tudatpy.estimation.observable_models_setup.model_settings.n_way_doppler_averaged_from_one_way_links

.. autofunction:: tudatpy.estimation.observable_models_setup.model_settings.dsn_n_way_doppler_averaged

.. autofunction:: tudatpy.estimation.observable_models_setup.model_settings.one_way_doppler_instantaneous

.. autofunction:: tudatpy.estimation.observable_models_setup.model_settings.two_way_doppler_instantaneous

.. autofunction:: tudatpy.estimation.observable_models_setup.model_settings.two_way_doppler_instantaneous_from_one_way_links

.. autofunction:: tudatpy.estimation.observable_models_setup.model_settings.angular_position

.. autofunction:: tudatpy.estimation.observable_models_setup.model_settings.relative_angular_position

.. autofunction:: tudatpy.estimation.observable_models_setup.model_settings.cartesian_position

.. autofunction:: tudatpy.estimation.observable_models_setup.model_settings.relative_cartesian_position

.. autofunction:: tudatpy.estimation.observable_models_setup.model_settings.cartesian_velocity

.. autofunction:: tudatpy.estimation.observable_models_setup.model_settings.euler_angles_313


Enumerations
------------
.. currentmodule:: tudatpy.estimation.observable_models_setup.model_settings

.. autosummary::

   ObservableType


.. autoclass:: tudatpy.estimation.observable_models_setup.model_settings.ObservableType
   :members:




Classes
-------
.. currentmodule:: tudatpy.estimation.observable_models_setup.model_settings

.. autosummary::

   DopplerProperTimeRateSettings

   ObservationModelSettings

   OneWayDopplerObservationModelSettings

   NWayRangeObservationModelSettings

.. autoclass:: tudatpy.estimation.observable_models_setup.model_settings.DopplerProperTimeRateSettings
   :members:

.. autoclass:: tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings
   :members:

.. autoclass:: tudatpy.estimation.observable_models_setup.model_settings.OneWayDopplerObservationModelSettings
   :members:

.. autoclass:: tudatpy.estimation.observable_models_setup.model_settings.NWayRangeObservationModelSettings

   :members:


.. _biases:

``biases``
==========

This module contains a set of factory functions for setting up the
biases for observation models.

The functions in this module create objects of type :class:`~tudatpy.estimation.biases.ObservationBiasSettings`,
which define settings for a type of observation bias. The main interface with Tudat is that these objects are used
as input to the observation model functions created in the :ref:`model_settings` module.

For an ideal observation :math:`h(t)`, this bias models created through the settings in this module
modify this into a biased (true) observation :math:`\tilde{h}(t)`, with the mapping from ideal to true
observation defined by the specific factory functions below.

More details on the link with the observation model is given on the `user guide <https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-model-setup.html>`_ .


Functions
---------
.. currentmodule:: tudatpy.estimation.observable_models_setup.biases

.. autosummary::

   absolute_bias

   relative_bias

   time_drift_bias

   time_bias

   combined_bias

   arcwise_absolute_bias

   arcwise_absolute_bias_per_time

   arcwise_relative_bias

   arcwise_relative_bias_per_time

   arc_wise_time_drift_bias

   arc_wise_time_drift_bias_per_time

   arcwise_time_bias

   clock_induced_bias

   two_way_time_scale_range_bias

.. autofunction:: tudatpy.estimation.observable_models_setup.biases.absolute_bias

.. autofunction:: tudatpy.estimation.observable_models_setup.biases.relative_bias

.. autofunction:: tudatpy.estimation.observable_models_setup.biases.combined_bias

.. autofunction:: tudatpy.estimation.observable_models_setup.biases.arcwise_absolute_bias

.. autofunction:: tudatpy.estimation.observable_models_setup.biases.arcwise_absolute_bias_per_time

.. autofunction:: tudatpy.estimation.observable_models_setup.biases.arcwise_relative_bias

.. autofunction:: tudatpy.estimation.observable_models_setup.biases.arcwise_relative_bias_per_time

.. autofunction:: tudatpy.estimation.observable_models_setup.biases.clock_induced_bias

.. autofunction:: tudatpy.estimation.observable_models_setup.biases.time_bias

.. autofunction:: tudatpy.estimation.observable_models_setup.biases.arcwise_time_bias

.. autofunction:: tudatpy.estimation.observable_models_setup.biases.time_drift_bias

.. autofunction:: tudatpy.estimation.observable_models_setup.biases.arc_wise_time_drift_bias

.. autofunction:: tudatpy.estimation.observable_models_setup.biases.arc_wise_time_drift_bias_per_time

.. autofunction:: tudatpy.estimation.observable_models_setup.biases.two_way_time_scale_range_bias


Classes
-------
.. currentmodule:: tudatpy.estimation.observable_models_setup.biases

.. autosummary::

   ObservationBiasSettings

.. autoclass:: tudatpy.estimation.observable_models_setup.biases.ObservationBiasSettings
   :members:
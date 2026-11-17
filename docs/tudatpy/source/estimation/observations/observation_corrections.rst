.. _observation_corrections:

``observation_corrections``
===========================

This module contains functionality to compute several observation corrections.
The corrections are requested on the tracking data with
:class:`~tudatpy.data_input.tracking_data.AngularObservationCorrectionSettings`, and computed once when the
observation dataset is created with ``apply_corrections=True``, before it is used in the estimation. Consequently, a reasonably accurate a priori ephemeris of the observed
object is required to compute the corrections accurately (e.g. retrieved from JPL Horizons).


Functions
---------

.. autosummary::

   ~tudatpy.estimation.observations.observation_corrections.photocenter_correction.photocenter_correction_angular_observations

   ~tudatpy.estimation.observations.observation_corrections.light_deflection_correction.light_deflection_correction_angular_observations

.. autofunction:: tudatpy.estimation.observations.observation_corrections.photocenter_correction.photocenter_correction_angular_observations

.. autofunction:: tudatpy.estimation.observations.observation_corrections.light_deflection_correction.light_deflection_correction_angular_observations

.. _observation_corrections:

``observation_corrections``
===========================

This module contains functionality to compute and apply several observation corrections.
These corrections are computed once and applied on an existing set of observations before
handing them to the estimation. As such, a reasonably accurate a priori ephemeris of the
object should be known to accurately compute the corrections. These can be retrieved from
e.g. JPL Horizons.


Functions
---------

.. autosummary::

   ~tudatpy.estimation.observations.observation_corrections.photocenter_correction.photocenter_correction_angular_observations

   ~tudatpy.estimation.observations.observation_corrections.photocenter_correction.apply_photocenter_correction_to_observation_collection

   ~tudatpy.estimation.observations.observation_corrections.light_deflection_correction.light_deflection_correction_angular_observations

   ~tudatpy.estimation.observations.observation_corrections.light_deflection_correction.apply_light_deflection_correction_to_observation_collection

.. autofunction:: tudatpy.estimation.observations.observation_corrections.photocenter_correction.photocenter_correction_angular_observations

.. autofunction:: tudatpy.estimation.observations.observation_corrections.photocenter_correction.apply_photocenter_correction_to_observation_collection

.. autofunction:: tudatpy.estimation.observations.observation_corrections.light_deflection_correction.light_deflection_correction_angular_observations

.. autofunction:: tudatpy.estimation.observations.observation_corrections.light_deflection_correction.apply_light_deflection_correction_to_observation_collection

.. _tracking_data_radar_utilities:

``radar_utilities``
===================

.. automodule:: tudatpy.data_input.tracking_data.radar_utilities

This submodule defines the canonical radar table shared by the MPC and JPL
readers. It provides conversion, filtering and validation functions, together
with the supplementary transmitter-frequency data needed by Tudat's radar
observation models.

.. currentmodule:: tudatpy.data_input.tracking_data.radar_utilities

Conversion
----------

.. autofunction:: radar_data_from_raw
.. autofunction:: radar_data_from_table
.. autofunction:: radar_data_to_tracking_data

Filtering and validation
------------------------

.. autofunction:: empty_radar_table
.. autofunction:: validate_radar_data
.. autofunction:: filter_radar_data

Frequency bands
---------------

.. autofunction:: radar_frequency_band_string_from_hz

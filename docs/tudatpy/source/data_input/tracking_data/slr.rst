.. _slr:

``slr``
=======

.. currentmodule:: tudatpy.data_input.tracking_data.slr

The CRD reader retains station/target headers, session correction flags,
system-specific transmit wavelengths, normal-point statistics, measured weather,
range supplements and calibration records. It does not apply atmospheric,
centre-of-mass, amplitude or calibration corrections. Inspect the H4 flags before
adding any correction to an observation model; a recorded calibration must not
be applied twice.

``extract_normal_point_measurements`` and ``extract_full_rate_measurements``
return an ordered list of ``(Time, one_way_range_metres)`` pairs. The variants
ending in ``_from_passes`` concatenate multiple sessions in input order. These
are lists rather than maps: coincident event epochs and distinct points that
round to the same float are retained. The native parser keeps high-resolution
event times when constructing :class:`~tudatpy.astro.time_representation.Time`.
The remaining scalar metadata fields use ordinary Python floats.

Extraction accepts raw two-way ranges with a recognised station UTC reference
(H2 codes 3, 4 or 7). Code 4 is UTC(GPS), not GPS coordinate time. The event epoch
is unchanged: the normal point's event flag still determines whether it denotes
transmission, reflection or reception. One-way range is ``c * TOF / 2``; this
unit conversion does not apply a clock transfer or a centre-of-mass correction.
Other range types remain available as raw records for explicit handling.

.. autosummary::

   read_slr_data
   extract_normal_point_measurements
   extract_normal_point_measurements_from_passes
   extract_full_rate_measurements
   extract_full_rate_measurements_from_passes
   convert_crd_two_way_time_of_flight_to_slr_range
   get_station_wavelengths
   CrdPass
   CrdPassConfigurationData
   CrdPassData
   CrdSystemConfigurationRecord
   CrdNormalPointRecord
   CrdFullRateRecord
   CrdMeteoRecord
   CrdRangeSupplementRecord
   CrdCalibrationRecord

.. automodule:: tudatpy.data_input.tracking_data.slr
   :members:

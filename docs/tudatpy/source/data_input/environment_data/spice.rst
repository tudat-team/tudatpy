.. _environment_data_spice:
.. _spice:

``spice``
=========

.. automodule:: tudatpy.data_input.environment_data.spice

This submodule contains Tudat's direct SPICE interface. It provides utilities
for loading and clearing SPICE kernels, converting SPICE time representations,
retrieving ephemeris and frame data, and querying body properties from the
loaded kernel pool. The old :mod:`tudatpy.interface.spice` module is kept only
as a compatibility shim; new code should import these functions from
:mod:`tudatpy.data_input.environment_data.spice`.

.. currentmodule:: tudatpy.data_input.environment_data.spice

Time Conversions
----------------

.. autosummary::

   convert_julian_date_to_ephemeris_time
   convert_ephemeris_time_to_julian_date
   convert_date_string_to_ephemeris_time
   get_approximate_utc_from_tdb

.. autofunction:: convert_julian_date_to_ephemeris_time
.. autofunction:: convert_ephemeris_time_to_julian_date
.. autofunction:: convert_date_string_to_ephemeris_time
.. autofunction:: get_approximate_utc_from_tdb

Ephemerides And Frames
----------------------

.. autosummary::

   get_body_cartesian_position_at_epoch
   get_body_cartesian_state_at_epoch
   get_cartesian_state_from_tle_at_epoch
   compute_rotation_matrix_between_frames
   compute_state_rotation_matrix_between_frames
   compute_rotation_matrix_derivative_between_frames
   get_angular_velocity_vector_of_frame_in_original_frame

.. autofunction:: get_body_cartesian_position_at_epoch
.. autofunction:: get_body_cartesian_state_at_epoch
.. autofunction:: get_cartesian_state_from_tle_at_epoch
.. autofunction:: compute_rotation_matrix_between_frames
.. autofunction:: compute_state_rotation_matrix_between_frames
.. autofunction:: compute_rotation_matrix_derivative_between_frames
.. autofunction:: get_angular_velocity_vector_of_frame_in_original_frame

Body Properties
---------------

.. autosummary::

   get_body_properties
   get_body_gravitational_parameter
   get_average_radius
   convert_body_name_to_naif_id
   convert_naif_id_to_body_name
   check_body_property_in_kernel_pool

.. autofunction:: get_body_properties
.. autofunction:: get_body_gravitational_parameter
.. autofunction:: get_average_radius
.. autofunction:: convert_body_name_to_naif_id
.. autofunction:: convert_naif_id_to_body_name
.. autofunction:: check_body_property_in_kernel_pool

Kernel Management
-----------------

.. autosummary::

   load_standard_kernels
   load_standard_deprecated_kernels
   load_kernel
   clear_kernels
   get_total_count_of_kernels_loaded
   continue_after_errors
   suppress_error_output

.. autofunction:: load_standard_kernels
.. autofunction:: load_standard_deprecated_kernels
.. autofunction:: load_kernel
.. autofunction:: clear_kernels
.. autofunction:: get_total_count_of_kernels_loaded
.. autofunction:: continue_after_errors
.. autofunction:: suppress_error_output

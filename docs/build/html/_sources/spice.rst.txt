``spice``
=========
Interface to the SPICE package.

This module provides an interface to NAIF's ``SPICE`` package.












Functions
---------
.. currentmodule:: tudatpy.interface.spice

.. autosummary::

   convert_julian_date_to_ephemeris_time

   convert_ephemeris_time_to_julian_date

   convert_date_string_to_ephemeris_time

   get_body_cartesian_state_at_epoch

   get_body_cartesian_position_at_epoch

   get_cartesian_state_from_tle_at_epoch

   compute_rotation_matrix_derivative_between_frames

   get_angular_velocity_vector_of_frame_in_original_frame

   get_body_properties

   get_body_gravitational_parameter

   get_average_radius

   convert_body_name_to_naif_id

   check_body_property_in_kernel_pool

   load_standard_kernels

   get_total_count_of_kernels_loaded

   load_kernel

   clear_kernels



.. autofunction:: tudatpy.interface.spice.convert_julian_date_to_ephemeris_time

.. autofunction:: tudatpy.interface.spice.convert_ephemeris_time_to_julian_date

.. autofunction:: tudatpy.interface.spice.convert_date_string_to_ephemeris_time

.. autofunction:: tudatpy.interface.spice.get_body_cartesian_state_at_epoch

.. autofunction:: tudatpy.interface.spice.get_body_cartesian_position_at_epoch

.. autofunction:: tudatpy.interface.spice.get_cartesian_state_from_tle_at_epoch

.. autofunction:: tudatpy.interface.spice.compute_rotation_matrix_derivative_between_frames

.. autofunction:: tudatpy.interface.spice.get_angular_velocity_vector_of_frame_in_original_frame

.. autofunction:: tudatpy.interface.spice.get_body_properties

.. autofunction:: tudatpy.interface.spice.get_body_gravitational_parameter

.. autofunction:: tudatpy.interface.spice.get_average_radius

.. autofunction:: tudatpy.interface.spice.convert_body_name_to_naif_id

.. autofunction:: tudatpy.interface.spice.check_body_property_in_kernel_pool

.. autofunction:: tudatpy.interface.spice.load_standard_kernels

.. autofunction:: tudatpy.interface.spice.get_total_count_of_kernels_loaded

.. autofunction:: tudatpy.interface.spice.load_kernel

.. autofunction:: tudatpy.interface.spice.clear_kernels






Classes
-------
.. currentmodule:: tudatpy.interface.spice

.. autosummary::

   SpiceEphemeris



.. autoclass:: tudatpy.interface.spice.SpiceEphemeris
   :members:
   :special-members: __init__






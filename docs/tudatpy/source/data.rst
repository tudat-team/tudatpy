.. _data:

``data``
========
Interfacing of Tudat(py) to and from other applications.

This module provides different functionalities that allow to export results from Tudat(py)
to other softwares such as MATLAB, to post-process results. A set of methods are also provided
to read inputs from other softwares and integrate them with Tudat(py).


.. toctree::
   :maxdepth: 2
   :caption: Modules

   data/horizons
   data/gaia
   data/mpc
   data/sbdb
   data/discos
   data/spacetrack
   data/processTrk234
   data/mission_data_downloader
   data/coma_model


Functions
---------
.. currentmodule:: tudatpy.data

.. autosummary::

   save2txt
   save_time_history_to_file
   read_vector_history_from_file
   read_matrix_history_from_file
   get_resource_path
   get_ephemeris_path
   get_earth_orientation_path
   get_quadrature_path
   get_spice_kernel_path
   get_atmosphere_tables_path
   get_gravity_models_path
   get_space_weather_path


   grail_antenna_file_reader

   grail_mass_level_0_file_reader

   grail_mass_level_1_file_reader

   read_fdets_file

   read_ifms_file

   read_odf_file

   read_solar_activity_data

   set_dsn_weather_data_in_ground_stations

   set_estrack_weather_data_in_ground_stations
.. autofunction:: tudatpy.data.save2txt

.. autofunction:: tudatpy.data.save_time_history_to_file

.. autofunction:: tudatpy.data.read_vector_history_from_file

.. autofunction:: tudatpy.data.read_matrix_history_from_file

.. autofunction:: tudatpy.data.get_resource_path

.. autofunction:: tudatpy.data.get_ephemeris_path

.. autofunction:: tudatpy.data.get_earth_orientation_path

.. autofunction:: tudatpy.data.get_quadrature_path

.. autofunction:: tudatpy.data.get_spice_kernel_path

.. autofunction:: tudatpy.data.get_atmosphere_tables_path

.. autofunction:: tudatpy.data.get_gravity_models_path

.. autofunction:: tudatpy.data.get_space_weather_path

.. autofunction:: tudatpy.data.grail_antenna_file_reader

.. autofunction:: tudatpy.data.grail_mass_level_0_file_reader

.. autofunction:: tudatpy.data.grail_mass_level_1_file_reader

.. autofunction:: tudatpy.data.read_fdets_file

.. autofunction:: tudatpy.data.read_ifms_file

.. autofunction:: tudatpy.data.read_odf_file

.. autofunction:: tudatpy.data.read_solar_activity_data

.. autofunction:: tudatpy.data.set_dsn_weather_data_in_ground_stations

.. autofunction:: tudatpy.data.set_estrack_weather_data_in_ground_stations

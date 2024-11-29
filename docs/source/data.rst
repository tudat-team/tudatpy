``data``
======
Interfacing of Tudat(py) to and from other applications.


This module provides with different functionalities that allow to export results from Tudat(py)
to other softwares such as MATLAB, to post-process results. A set of methods are also provided
to read inputs from other softwares and integrate them with Tudat(py).


.. toctree::
   :maxdepth: 2
   :caption: Modules

   mpc




Functions
---------
.. currentmodule:: tudatpy.data

.. autosummary::

   save2txt

   save_time_history_to_file

   get_resource_path

   get_ephemeris_path

   get_earth_orientation_path

   get_quadrature_path

   get_spice_kernel_path

   get_atmosphere_tables_path

   get_gravity_models_path

   get_space_weather_path

   read_vector_history_from_file

   read_matrix_history_from_file



.. autofunction:: tudatpy.data.save2txt

.. autofunction:: tudatpy.data.save_time_history_to_file

.. autofunction:: tudatpy.data.get_resource_path

.. autofunction:: tudatpy.data.get_ephemeris_path

.. autofunction:: tudatpy.data.get_earth_orientation_path

.. autofunction:: tudatpy.data.get_quadrature_path

.. autofunction:: tudatpy.data.get_spice_kernel_path

.. autofunction:: tudatpy.data.get_atmosphere_tables_path

.. autofunction:: tudatpy.data.get_gravity_models_path

.. autofunction:: tudatpy.data.get_space_weather_path

.. autofunction:: tudatpy.data.read_vector_history_from_file

.. autofunction:: tudatpy.data.read_matrix_history_from_file









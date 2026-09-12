.. _data_input:

``data_input``
==============

This submodule contains functionality for retrieving and reading data from
external sources. This includes user-provided files, Tudat resource files,
mission archives, and external services or packages.

The functionality in this module is limited to getting data into Python and
representing it in a form that Tudat can use. The use of these data in an
environment, propagation, or estimation model is handled by the corresponding
modules in :ref:`dynamics` and :ref:`estimation`.


* :ref:`environment_data`:  Interfaces to various file formats/packages/services (includes Spice, JPL Horizons, Small-Body Database and Spacetrack) that provide data that can be used in Tudat to define physical properties of solar system bodies, both natural and artificial.  These interfaces are often used by functions in
  :ref:`environment_setup` to define settings for Tudat body objects.
* :ref:`tracking_data`: Interfaces to various file formats/packages/services that contain tracking data
  that can be used in orbit and parameter estimation, including radio tracking data from DSN and ESTRACK and optical
  data from the MPC. This module loads tracking data and associated metadata; conversion to Tudat-compatible data
  containers used in estimation is handled in the :ref:`observations` module.
* :ref:`data_retrieval`: Download utilities for
  external data products, including mission-archive files, SPICE kernels,
  radio-science files, and media-correction data such as IONEX and
  VMF files.
* :ref:`resource_paths`: Functions that return
  paths to Tudat resource folders and common input-data locations, such as
  ephemerides, gravity models, atmosphere tables, Earth-orientation data, and
  test data.

.. toctree::
   :maxdepth: 2
   :caption: Modules

   data_input/environment_data
   data_input/tracking_data
   data_input/data_retrieval
   data_input/resource_paths


.. automodule:: tudatpy.data_input
   :members:

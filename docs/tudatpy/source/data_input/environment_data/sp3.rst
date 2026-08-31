.. _environment_data_sp3:

``sp3``
========

This submodule reads SP3 precise-orbit files and exposes their states and
metadata. It does not perform reference-frame or time-scale conversions.
Use the SP3 factories in
:mod:`tudatpy.dynamics.environment_setup.ephemeris` to create ephemeris
settings and apply supported frame transformations.

.. automodule:: tudatpy.data_input.environment_data.sp3

.. currentmodule:: tudatpy.data_input.environment_data.sp3

Functions
---------

.. autosummary::

   read_sp3_file

.. autofunction:: read_sp3_file

Classes
-------

.. autosummary::

   Sp3FileContents

.. autoclass:: Sp3FileContents
   :members:

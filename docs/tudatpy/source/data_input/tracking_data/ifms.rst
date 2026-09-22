.. _tracking_data_ifms:

``ifms``
========

.. automodule:: tudatpy.data_input.tracking_data.ifms

This submodule contains functionality to load tracking data from IFMS files.
IFMS files contain radio tracking data, typically closed-loop Doppler data, for
a number of ESA deep space missions, optionally with station tropospheric
corrections that can be applied during loading. The IFMS-to-OCC interface and
file format are described in :cite:t:`ifmsOccFtp2006`. The
:func:`read_ifms_data` function is the main interface for loading the data and
converting it to objects that Tudat can process further; see also
:ref:`tracking_data`. All other functionality in this module is
reserved for better understanding what data is being loaded, and in some cases
manipulating it, before it is processed into Tudat-compatible objects.

.. currentmodule:: tudatpy.data_input.tracking_data.ifms

.. autofunction:: read_ifms_data

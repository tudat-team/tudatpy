.. _tracking_data_fdets:

``fdets``
=========

.. automodule:: tudatpy.data_input.tracking_data.fdets

This submodule contains functionality to load tracking data from Fdets files.
Fdets files contain Doppler frequency observables and associated metadata such
as signal-to-noise ratio, spectral maximum, and Doppler noise for observations
made using the PRIDE experiment. The :func:`read_fdets_data` function is the
main interface for loading the data and converting it to objects that Tudat can
process further; see also
:ref:`tracking_data`. All other functionality in this module is
reserved for better understanding what data is being loaded, and in some cases
manipulating it, before it is processed into Tudat-compatible objects.

.. currentmodule:: tudatpy.data_input.tracking_data.fdets

.. autofunction:: read_fdets_data

Supporting API
--------------

The class below defines how dates are represented in the input Fdets files. It
is used when calling :func:`read_fdets_data`; no separate parsing or conversion
step is normally needed.

.. autoclass:: FdetDateFormat

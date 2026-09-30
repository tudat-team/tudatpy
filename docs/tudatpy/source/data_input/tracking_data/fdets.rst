.. _tracking_data_fdets:

``fdets``
=========

.. automodule:: tudatpy.data_input.tracking_data.fdets

This submodule contains functionality to load tracking data from FDETS files, as produced by the PRIDE experiment (see the `data processing description <https://doi.org/10.1017/pasa.2021.56>`_).
FDETS files contain open-loop Doppler frequency observables and associated metadata such
as signal-to-noise ratio, spectral maximum, and Doppler noise for observations
made using the PRIDE experiment.
To use these data in an orbit estimation, the transmitter frequency needs to be defined, or loaded from another data source (such as ODF or IFMS files).

The :func:`read_fdets_data` function is the
main interface for loading the data and converting it to objects that Tudat can
process further; see also
:ref:`tracking_data`. All other functionality in this module is
reserved for better understanding what data is being loaded, and in some cases
manipulating it, before it is processed into Tudat-compatible objects.

.. currentmodule:: tudatpy.data_input.tracking_data.fdets

.. autofunction:: read_fdets_data

Supporting API
--------------

The class below defines how dates are represented in the input FDETS files. It
is used when calling :func:`read_fdets_data`; no separate parsing or conversion
step is normally needed.

.. autoclass:: FdetDateFormat

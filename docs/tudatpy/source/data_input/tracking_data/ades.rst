.. _tracking_data_ades:

``ades``
===============

.. automodule:: tudatpy.data_input.tracking_data.ades

This submodule loads optical astrometry, space-based astrometry and radar data
from MPC ADES files, in both psv and xml formats. The supported record layouts are
described by the Minor Planet Center for
`optical, radar, stellar occultation and offset observations <https://minorplanetcenter.net/iau/info/ADES.html>`_. The
:func:`read_ades_data` function is the main interface for loading the data
and converting it to objects that Tudat can process further; see also
:ref:`tracking_data`. Offset observations, radar observations that are not
w.r.t. the centre of mass, and observations from Earth-based roving observers
are not yet supported, and are skipped with a warning.

.. currentmodule:: tudatpy.data_input.tracking_data.ades

.. autofunction:: read_ades_data

Supporting API
--------------

The function below is used internally by :func:`read_ades_data` to parse the ADES
records into an astropy table. It does not create Tudat
tracking-data objects, but can be used to inspect parsed optical, radar, stellar occultations and offset data.
Radar observations are available in the returned table's metadata as a canonical radar table.
Ability to read stellar occultations and offset data is currently being implemented.

.. autofunction:: parse_ades_file

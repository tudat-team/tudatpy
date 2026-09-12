.. _tracking_data_mpc:

``mpc``
=======

.. automodule:: tudatpy.data_input.tracking_data.mpc

This submodule contains functionality to retrieve optical tracking data from
the Minor Planet Center (MPC) database for asteroids and comets. The
:class:`BatchMPC` class wraps the MPC interface provided by ``astroquery`` and
adds Tudat-specific processing, including optional observation weights based on
:cite:t:`veres2017` and star-catalog bias corrections based on
:cite:t:`eggl2020`. The :func:`read_mpc_data` function is the main interface
for loading the data and converting it to objects that Tudat can process
further; see also
:ref:`tracking_data`. All other functionality in this module is
reserved for better understanding what data is being loaded, and in some cases
manipulating it, before it is processed into Tudat-compatible objects.

.. currentmodule:: tudatpy.data_input.tracking_data.mpc

.. autofunction:: read_mpc_data

Supporting API
--------------

The class below exposes the intermediate MPC batch used by
:func:`read_mpc_data`. It can be used directly to retrieve observations,
inspect the resulting table and metadata, and then convert the batch to Tudat
tracking-data objects.

.. autoclass:: BatchMPC
   :members:

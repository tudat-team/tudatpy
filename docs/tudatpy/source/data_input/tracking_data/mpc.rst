.. _tracking_data_mpc:

``mpc``
=======

.. automodule:: tudatpy.data_input.tracking_data.mpc

This submodule retrieves optical, space-based and radar tracking data from the
Minor Planet Center (MPC) database for asteroids and comets. The
:class:`BatchMPC` class wraps the MPC interface provided by ``astroquery`` and
adds Tudat-specific processing, including optional optical observation weights
based on :cite:t:`veres2017` and star-catalog bias corrections based on
:cite:t:`eggl2020`. Raw MPC 80-column retrieval is available through the
``use_mpc80_format`` argument and is required for space-based and radar
observations. The :func:`read_mpc_data` function is the main interface for
loading the data and converting it to objects that Tudat can process further;
see also :ref:`tracking_data`.

.. currentmodule:: tudatpy.data_input.tracking_data.mpc

.. autofunction:: read_mpc_data

Supporting API
--------------

The class below exposes the intermediate MPC batch used by
:func:`read_mpc_data`. It can be used directly to retrieve observations,
inspect the resulting table and metadata, and then convert the batch to Tudat
tracking-data objects. It also provides a method for converting spacecraft
positions from space-based observations into state histories.

.. autoclass:: BatchMPC
   :members:

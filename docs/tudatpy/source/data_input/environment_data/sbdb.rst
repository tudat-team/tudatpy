.. _environment_data_sbdb:

``sbdb``
========

This module contains a wrapper for selected data from JPL's Small-Body
Database (SBDB). The current Tudat interface uses only a limited part of the
available SBDB information: object names, SPK identifiers, diameters, and
gravitational parameters when available.

For direct environment setup, SBDB data are used by
:func:`~tudatpy.dynamics.environment_setup.gravity_field.sbdb_wrapper.central_sbdb` and
:func:`~tudatpy.dynamics.environment_setup.gravity_field.sbdb_wrapper.central_sbdb_density`
to create point-mass gravity field settings for small bodies. The
:class:`~tudatpy.data_input.environment_data.sbdb.SBDBquery` class can be used directly when these physical properties
need to be inspected before constructing the environment.

.. automodule:: tudatpy.data_input.environment_data.sbdb

Classes
-------
.. currentmodule:: tudatpy.data_input.environment_data.sbdb

.. autosummary::

   SBDBquery

.. autoclass:: SBDBquery
   :members:

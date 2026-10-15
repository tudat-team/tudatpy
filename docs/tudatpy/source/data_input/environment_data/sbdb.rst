.. _environment_data_sbdb:

``sbdb``
========

This module provides two independent interfaces to JPL's Small-Body Database
(SBDB). :class:`~tudatpy.data_input.environment_data.sbdb.SBDBquery` queries one
object and exposes its properties.
:class:`~tudatpy.data_input.environment_data.sbdb.SBDBbatch` downloads or loads
the catalogue as a table; its ``get`` method filters that table by primary MPC
designation without making another request.
The interface is based on the `astroquery <https://github.com/astropy/astroquery>`_ :cite:`ginsburg2019astroquery` Python package.


For example, a catalogue query can select an object and save the table for
later reuse:

.. code-block:: python

    from tudatpy.data_input.environment_data.sbdb import SBDBbatch

    batch = SBDBbatch(fields=["pdes", "name"])
    vesta = batch.get(4)
    batch.to_csv("sbdb.csv")
    reloaded = SBDBbatch(csv_path="sbdb.csv")

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

   SBDBbatch
   SBDBquery

.. autoclass:: SBDBbatch
   :members:

.. autoclass:: SBDBquery
   :members:

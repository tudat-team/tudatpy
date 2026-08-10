.. _environment_data_horizons:

``horizons``
============

This module contains an interface to JPL Horizons for retrieving ephemeris data
of natural and artificial Solar System bodies. In Tudat, these data are most
often used to define ephemeris settings through
:func:`~tudatpy.dynamics.environment_setup.ephemeris.horizons_wrapper.jpl_horizons_from_query`
or :func:`~tudatpy.dynamics.environment_setup.ephemeris.horizons_wrapper.jpl_horizons`.
:func:`~tudatpy.dynamics.environment_setup.ephemeris.horizons_wrapper.jpl_horizons_from_query`
is equivalent to :func:`~tudatpy.dynamics.environment_setup.ephemeris.horizons_wrapper.jpl_horizons`,
but allows the :class:`~tudatpy.data_input.environment_data.horizons.HorizonsQuery`
to be inspected or modified before conversion to ephemeris settings.

The :class:`~tudatpy.data_input.environment_data.horizons.HorizonsQuery` class can also be used directly to query Horizons
and inspect the returned vectors. The :class:`~tudatpy.data_input.environment_data.horizons.HorizonsBatch` class is a
convenience wrapper for applying the same query settings to multiple bodies.

.. automodule:: tudatpy.data_input.environment_data.horizons

Classes
-------
.. currentmodule:: tudatpy.data_input.environment_data.horizons

.. autosummary::

   HorizonsQuery
   HorizonsBatch

.. autoclass:: HorizonsQuery
   :members:

.. autoclass:: HorizonsBatch
   :members:

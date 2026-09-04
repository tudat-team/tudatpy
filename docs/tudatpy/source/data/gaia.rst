.. _gaia:

``gaia``
========
This module contains interfaces to retrieve and process Gaia solar-system astrometry for asteroids and
Gaia-derived asteroid orbits and covariance data from the Gaia archives.


Classes
-------
.. currentmodule:: tudatpy.data.gaia

.. autosummary::

   GaiaAstrometry

.. autoclass:: tudatpy.data.gaia.GaiaAstrometry
   :members:
   :special-members: __init__

Functions
---------
.. autosummary::

   get_state_from_gaia_archive

   get_state_covariance_from_gaia_archive

   get_kepler_covariance_from_gaia_archive

   gaia_object_catalog

   generate_astrometry_parquet

   generate_asteroid_parquet


.. autofunction:: tudatpy.data.gaia.get_state_from_gaia_archive

.. autofunction:: tudatpy.data.gaia.get_state_covariance_from_gaia_archive

.. autofunction:: tudatpy.data.gaia.get_kepler_covariance_from_gaia_archive

.. autofunction:: tudatpy.data.gaia.gaia_object_catalog

.. autofunction:: tudatpy.data.gaia.generate_astrometry_parquet

.. autofunction:: tudatpy.data.gaia.generate_asteroid_parquet

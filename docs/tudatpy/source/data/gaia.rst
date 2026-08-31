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

   GaiaAsteroids


.. autoclass:: tudatpy.data.gaia.GaiaAstrometry
   :members:
   :special-members: __init__

.. autoclass:: tudatpy.data.gaia.GaiaAsteroids
   :members:
   :special-members: __init__


Functions
---------
.. autosummary::

   generate_astrometry_parquet

   generate_asteroid_parquet


.. autofunction:: tudatpy.data.gaia.generate_astrometry_parquet

.. autofunction:: tudatpy.data.gaia.generate_asteroid_parquet

.. _environment_data_damit:

``damit``
=========

This module contains an interface to the Database of Asteroid Models from
Inversion Techniques (DAMIT). This database contains polyhedron shape models
and rotational parameters for many asteroids. Through this interface, these
models can be queried and automatically downloaded, and then converted into
shape models and rotational ephemeris settings. It can also approximate a best
fitting ellipsoid on the shape models.

See https://damit.cuni.cz/projects/damit/pages/documentation for documentation.

.. automodule:: tudatpy.data_input.environment_data.damit

Functions
---------
.. currentmodule:: tudatpy.data_input.environment_data.damit

.. autosummary::

   damit_model_query

.. autofunction:: damit_model_query

Classes
-------
.. currentmodule:: tudatpy.data_input.environment_data.damit

.. autosummary::

   DamitModel

.. autoclass:: DamitModel
   :members:

.. _environment_data_damit:

``damit``
=========

This module contains an interface to the Database of Asteroid Models from
Inversion Techniques (DAMIT). This database contains polyhedron shape models
and rotational parameters for many asteroids. Through this interface, these
models can be queried and converted into rotational models and shape models.
When attempting to use the module for the first time, you will be asked to download
the DAMIT catalog files. This is a one time download of approximately 4 GB.

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

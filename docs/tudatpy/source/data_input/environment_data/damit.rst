.. _environment_data_damit:

``damit``
=========

This module contains an interface to the Database of Asteroid Models from
Inversion Techniques (DAMIT). This database contains polyhedron shape models
and rotational parameters for many asteroids. Through this interface, these
models can be queried and automatically downloaded. It can also approximate a
best fitting ellipsoid on the shape models. This module does not support
tumblers at the moment.

To convert a model into shape and rotation model settings, use
:func:`~tudatpy.dynamics.environment_setup.shape.damit_wrapper.damit_polyhedron`
and :func:`~tudatpy.dynamics.environment_setup.rotation_model.damit_wrapper.damit_polyhedron`
(or the ``damit_best_fitting_ellipsoid`` versions of both).

See https://damit.cuni.cz/projects/damit/pages/documentation for documentation.

.. automodule:: tudatpy.data_input.environment_data.damit

Functions
---------
.. currentmodule:: tudatpy.data_input.environment_data.damit

.. autosummary::

   list_damit_models

.. autofunction:: list_damit_models

Classes
-------
.. currentmodule:: tudatpy.data_input.environment_data.damit

.. autosummary::

   DamitModel

.. autoclass:: DamitModel
   :members:

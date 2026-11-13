.. _environment_data_damit:

``damit``
=========

This module provides an interface to the Database of Asteroid Models from
Inversion Techniques (`DAMIT <https://damit.cuni.cz/projects/damit/>`_)
:cite:p:`Durech2010`, which contains polyhedron shape models and rotation
states of asteroids.
:class:`~tudatpy.data_input.environment_data.damit.DamitModel` loads a model
from the database by MPC number or model ID (downloading it on first use), or
from files, and computes the tri-axial ellipsoid with the same moments of
inertia. Models of tumbling asteroids are not supported.

For example, the shape and rotation of Eros are set up from its DAMIT model by:

.. code-block:: python

    from tudatpy.data_input.environment_data.damit import DamitModel

    eros = DamitModel.query(mpc_number=433)
    body_settings.get("Eros").shape_settings = environment_setup.shape.damit_polyhedron(eros)
    body_settings.get("Eros").rotation_model_settings = environment_setup.rotation_model.damit(
        eros, "ECLIPJ2000", "Eros_fixed"
    )

:func:`~tudatpy.dynamics.environment_setup.shape.damit_best_fitting_ellipsoid`
gives the best fitting ellipsoid instead of the polyhedron, in the same
body-fixed frame. The DAMIT conventions are described in its
`documentation <https://damit.cuni.cz/projects/damit/pages/documentation>`_.

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

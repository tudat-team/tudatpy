.. _environment_data_discos:

``discos``
==========

This module contains an interface to ESA's DISCOS system for retrieving
properties of resident space objects. The retrieved properties can be used to
define body properties in :ref:`environment_setup`, for instance constant
aerodynamic coefficients through :func:`~tudatpy.dynamics.environment_setup.aerodynamic_coefficients.constant`
or mass properties through :func:`~tudatpy.dynamics.environment_setup.rigid_body.custom_mass_dependent_rigid_body_properties`.

To use this interface, you must have access to ESA's `DISCOSweb <https://discosweb.esoc.esa.int/>`_ service.

.. automodule:: tudatpy.data_input.environment_data.discos

Classes
-------
.. currentmodule:: tudatpy.data_input.environment_data.discos

.. autosummary::

   DiscosQuery

.. autoclass:: DiscosQuery
   :members:

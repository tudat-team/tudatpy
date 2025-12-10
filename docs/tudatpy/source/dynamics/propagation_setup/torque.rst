.. _torque:

``torque``
==========
Here, you will find a list of all torque models available in Tudat. The functions here all create settings
objects for torques that are provided as input to the :func:`~tudatpy.dynamics.propagation_setup.create_torque_models`
function. For more details on how these torques interface with the rest of the Tudat propagation framework,
see our user guide page on `rotational dynamics <https://docs.tudat.space/en/latest/user-guide/state-propagation/propagation-setup/rotational.html>`_

The functions in this submodule create objects of type :class:`~TorqueSettings` or (in case more information than only
the type of torque is needed to create the torque), one of its derived classes.


Functions
---------
.. currentmodule:: tudatpy.dynamics.propagation_setup.torque

.. autosummary::

   aerodynamic

   spherical_harmonic_gravitational

   second_degree_gravitational

   custom_torque



.. autofunction:: tudatpy.dynamics.propagation_setup.torque.aerodynamic

.. autofunction:: tudatpy.dynamics.propagation_setup.torque.spherical_harmonic_gravitational

.. autofunction:: tudatpy.dynamics.propagation_setup.torque.second_degree_gravitational

.. autofunction:: tudatpy.dynamics.propagation_setup.torque.custom_torque




Enumerations
------------
.. currentmodule:: tudatpy.dynamics.propagation_setup.torque

.. autosummary::

   AvailableTorque



.. autoclass:: tudatpy.dynamics.propagation_setup.torque.AvailableTorque
   :members:




Classes
-------
.. currentmodule:: tudatpy.dynamics.propagation_setup.torque

.. autosummary::

   TorqueSettings

   SphericalHarmonicTorqueSettings



.. autoclass:: tudatpy.dynamics.propagation_setup.torque.TorqueSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.torque.SphericalHarmonicTorqueSettings
   :members:




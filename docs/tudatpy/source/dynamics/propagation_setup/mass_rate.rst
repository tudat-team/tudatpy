.. _mass_rate:

``mass_rate``
=============
Here, you will find a list of all mass rate models available in Tudat, which are needed when propagating the mass of
a body numerically (for instance for a spacecraft under thrust). Some more specifics on this, and the overall interaction with Tudat, can be found on
the `user guide <https://docs.tudat.space/en/latest/user-guide/state-propagation/propagation-setup/mass.html>`_ .

The functions here all create settings
objects for mass rate models that are provided as input to the :func:`~tudatpy.dynamics.propagation_setup.create_mass_rate_models`
function. The functions in this submodule create objects of type :class:`~MassRateModelSettings` or (in case more information than only
the type of mass rate is needed to create the model), one of its derived classes.

Functions
---------
.. currentmodule:: tudatpy.dynamics.propagation_setup.mass_rate

.. autosummary::

   from_thrust

   custom_mass_rate



.. autofunction:: tudatpy.dynamics.propagation_setup.mass_rate.from_thrust

.. autofunction:: tudatpy.dynamics.propagation_setup.mass_rate.custom_mass_rate




Enumerations
------------
.. currentmodule:: tudatpy.dynamics.propagation_setup.mass_rate

.. autosummary::

   AvailableMassRateModels



.. autoclass:: tudatpy.dynamics.propagation_setup.mass_rate.AvailableMassRateModels
   :members:




Classes
-------
.. currentmodule:: tudatpy.dynamics.propagation_setup.mass_rate

.. autosummary::

   MassRateModelSettings

   FromThrustMassRateSettings

   CustomMassRateSettings



.. autoclass:: tudatpy.dynamics.propagation_setup.mass_rate.MassRateModelSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.mass_rate.FromThrustMassRateSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.mass_rate.CustomMassRateSettings
   :members:




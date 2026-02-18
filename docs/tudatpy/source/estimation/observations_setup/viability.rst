.. _viability:

``viability``
=============


In many cases, whether an observation at a given time should be realized will depend on a number of constraints that must be satisfied. We have termed such constraints 'observation viability conditions' in Tudat. The functions in this module are used create objects of type :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` that contain settings for these viability conditions. The main interface with Tudat is that these objects are can be added to :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects, which are in turn used as input to observation simulation (see :ref:`observations_simulation_settings`). The ancillary settings can be set in these objects directly using the :attr:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings.viability_settings_list` attribute, or by using the :func:`~tudatpy.estimation.observations_setup.viability.add_viability_check_to_all`, :func:`~tudatpy.estimation.observations_setup.viability.add_viability_check_to_observable` and :func:`~tudatpy.estimation.observations_setup.viability.add_viability_check_to_observable_for_link_ends` functions (which are functions of convenience to add the same ancillary settings to a set of observation simulation settings)

As an example, an ``observation_simulation_settings_list`` can be modified such that only observations above a 15 degree elevation angle at New Norcia (NNO) are accepted:

.. code-block:: python

    from tudatpy.estimation.observations_setup import viability
    from tudatpy.estimation.observable_models_setup import links
    import numpy as np

    # Create list of observation simulation settings (list[observations_simulation_settings.ObservationSimulationSettings])
    observation_simulation_settings_list = ...

    station_id = links.body_reference_point_link_end_id("Earth", "NNO")
    viability_settings_list = list()
    viability_settings_list.append(
        viability.elevation_angle_viability(
            station_id,
            np.deg2rad(15.0)
        )
    )
    viability.add_viability_check_to_all(
        observation_simulation_settings_list,
        viability_settings_list
    )

Functions
---------
.. currentmodule:: tudatpy.estimation.observations_setup.viability

.. autosummary::

   elevation_angle_viability

   body_avoidance_viability

   body_occultation_viability

   elevation_angle_viability_list

   body_avoidance_viability_list

   body_occultation_viability_list

   add_viability_check_to_all

   add_viability_check_to_observable

   add_viability_check_to_observable_for_link_ends

.. autofunction:: tudatpy.estimation.observations_setup.viability.elevation_angle_viability

.. autofunction:: tudatpy.estimation.observations_setup.viability.body_avoidance_viability

.. autofunction:: tudatpy.estimation.observations_setup.viability.body_occultation_viability

.. autofunction:: tudatpy.estimation.observations_setup.viability.elevation_angle_viability_list

.. autofunction:: tudatpy.estimation.observations_setup.viability.body_avoidance_viability_list

.. autofunction:: tudatpy.estimation.observations_setup.viability.body_occultation_viability_list

.. autofunction:: tudatpy.estimation.observations_setup.viability.add_viability_check_to_all

.. autofunction:: tudatpy.estimation.observations_setup.viability.add_viability_check_to_observable

.. autofunction:: tudatpy.estimation.observations_setup.viability.add_viability_check_to_observable_for_link_ends

Enumerations
------------
.. currentmodule:: tudatpy.estimation.observations_setup.viability

.. autosummary::

   ObservationViabilityType


.. autoclass:: tudatpy.estimation.observations_setup.viability.ObservationViabilityType
   :members:




Classes
-------
.. currentmodule:: tudatpy.estimation.observations_setup.viability

.. autosummary::

   ObservationViabilitySettings

.. autoclass:: tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings
   :members:
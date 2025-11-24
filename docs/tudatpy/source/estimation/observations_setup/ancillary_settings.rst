.. _ancillary_settings:

``ancillary_settings``
=======================

Some observation models require or can use additional quantities in addition to the time tag as input to the calculation of an observable, such a retransmission delay or a frequency band. In Tudat, we call these quantities ancillary settings. The functions in this module are used create objects of type :class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings` that contain these quantities. The full list of types on ancillary data that can be stored in this object is defined by the :class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationVariable` enum. The main interface with Tudat is that these objects are can be added to :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects, which are in turn used as input to observation simulation (see :ref:`observations_simulation_settings`). The ancillary settings can be set in these objects directly using the :attr:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings.ancillary_settings` attribute, or by using the :func:`~tudatpy.estimation.observations_setup.ancillary_settings.add_ancilliary_settings_to_observable` and :func:`~tudatpy.estimation.observations_setup.ancillary_settings.add_ancilliary_settings_to_observable_for_link_ends` functions (which are functions of convenience to add the same ancillary settings to a set of observation simulation settings)

This module provides two manners in which to create ancillary settings objects:

* A list of functions for specific observation types (:func:`~tudatpy.estimation.observations_setup.ancillary_settings.two_way_range_ancilliary_settings`, :func:`~tudatpy.estimation.observations_setup.ancillary_settings.two_way_doppler_ancilliary_settings`, etc.) that automatically creates an ancillary settings object with the type of settings relevant for a given observable type. In most cases, this is the easiest manner in which to define settings.

.. code-block:: python

     # Code snippet to show the creation of an ObservationAncillarySimulationSettings object
     from tudatpy.estimation.observations_setup import ancillary_settings

     # Example 1: Create ObservationAncillarySimulationSettings object using ancillary_settings.n_way_range_ancilliary_settings function
     # In this case the frequency bands of the retransmitter - we set it to x band.
     n_way_range_ancillary_settings = ancillary_settings.n_way_range_ancilliary_settings(frequency_bands=[ancillary_settings.FrequencyBands.x_band])

     # Show that this is indeed an ObservationAncillarySimulationSettings object
     print(n_way_range_ancillary_settings)

     # Example 2: Create ObservationAncillarySimulationSettings object using ancillary_settings.doppler_ancilliary_settings function
     # In this case the integration time (in seconds) has to be given as input - we set it to 60s
     doppler_ancillary_settings = ancillary_settings.doppler_ancilliary_settings(60)

     # Show that this is indeed an ObservationAncillarySimulationSettings object
     print(doppler_ancillary_settings)

* Manual creation of an empty :class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings` object, followed by the manual setting of all the relevant quantities using the :meth:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings.set_float_settings` and :meth:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings.set_float_list_settings` methods

Functions
---------
.. currentmodule:: tudatpy.estimation.observations_setup.ancillary_settings

.. autosummary::

   doppler_ancilliary_settings

   two_way_range_ancilliary_settings

   two_way_doppler_ancilliary_settings

   n_way_range_ancilliary_settings

   n_way_doppler_ancilliary_settings

   dsn_n_way_range_ancilliary_settings

   dsn_n_way_doppler_ancilliary_settings

   doppler_measured_frequency_ancillary_settings

   add_ancilliary_settings_to_observable

   add_ancilliary_settings_to_observable_for_link_ends

   dsn_default_turnaround_ratios

   cassini_turnaround_ratios

.. autofunction:: tudatpy.estimation.observations_setup.ancillary_settings.doppler_ancilliary_settings

.. autofunction:: tudatpy.estimation.observations_setup.ancillary_settings.two_way_range_ancilliary_settings

.. autofunction:: tudatpy.estimation.observations_setup.ancillary_settings.two_way_doppler_ancilliary_settings

.. autofunction:: tudatpy.estimation.observations_setup.ancillary_settings.n_way_range_ancilliary_settings

.. autofunction:: tudatpy.estimation.observations_setup.ancillary_settings.n_way_doppler_ancilliary_settings

.. autofunction:: tudatpy.estimation.observations_setup.ancillary_settings.dsn_n_way_range_ancilliary_settings

.. autofunction:: tudatpy.estimation.observations_setup.ancillary_settings.dsn_n_way_doppler_ancilliary_settings

.. autofunction:: tudatpy.estimation.observations_setup.ancillary_settings.doppler_measured_frequency_ancillary_settings

.. autofunction:: tudatpy.estimation.observations_setup.ancillary_settings.add_ancilliary_settings_to_observable

.. autofunction:: tudatpy.estimation.observations_setup.ancillary_settings.add_ancilliary_settings_to_observable_for_link_ends

.. autofunction:: tudatpy.estimation.observations_setup.ancillary_settings.dsn_default_turnaround_ratios

.. autofunction:: tudatpy.estimation.observations_setup.ancillary_settings.cassini_turnaround_ratios



Enumerations
------------
.. currentmodule:: tudatpy.estimation.observations_setup.ancillary_settings

.. autosummary::

   ObservationAncilliarySimulationVariable

   ObservationIntermediateSimulationVariable

   FrequencyBands


.. autoclass:: tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationVariable
   :members:

.. autoclass:: tudatpy.estimation.observations_setup.ancillary_settings.ObservationIntermediateSimulationVariable
   :members:

.. autoclass:: tudatpy.estimation.observations_setup.ancillary_settings.FrequencyBands
   :members:



Classes
-------
.. currentmodule:: tudatpy.estimation.observations_setup.ancillary_settings

.. autosummary::

    ObservationAncilliarySimulationSettings


.. autoclass:: tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings
   :members:
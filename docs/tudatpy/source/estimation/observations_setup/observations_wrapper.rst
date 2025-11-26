.. _observations_wrapper:

``observations_wrapper``
========================

The functionality in this module is used to either load or simulate observations, typically for later use in an orbit estimation. The functionality here is for loading/simulating data using Tudat-native functionality. We also provide functionality to load data from various sources using external packages/libraries (see the :ref:`data` module, including for TNF files from the DSN, and for astrometric data from the MPC). More information on using real data for an estimation in Tudat is given on the `user guide <https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-simulation/creating-observations/loading-real-data.html>`_

This module supports the loading of the following data types:

* ODF files: Data files produced by NASA Deep Space Network (DSN) according to the  `TRK 2-18 format <https://pds-geosciences.wustl.edu/radiosciencedocs/urn-nasa-pds-radiosci_documentation/DSN_TRK-2-18/dsn_trk-2-18.2008-02-29.pdf>`_. ODF files are binary files, which can be read by Tudat using the :func:`~tudatpy.estimation.observations_setup.observations_wrapper.process_odf_data_multiple_files` and converted into an :class:`~tudatpy.estimation.observations.ObservationCollection` using the :func:`~tudatpy.estimation.observations_setup.observations_wrapper.create_odf_observed_observation_collection` function. The ODF files also contain information on station transmission frequencies (ramp tables) which can be automatically set in the ground station objects in Tudat from the data in the ODF files using the :func:`~tudatpy.estimation.observations_setup.observations_wrapper.set_odf_information_in_bodies` function. As a function of convenience, the :func:`~tudatpy.estimation.observations_setup.observations_from_odf_files` takes care of each of these steps. At the moment, only two- and three-way Doppler data, and two-way sequential ranging data are supported in the file reading of ODF files
* IFMS files: Data files produced by ESA's ESTRACK tracking network, see `generic documentation <https://archives.esac.esa.int/psa/ftp/MARS-EXPRESS/MRS/MEX-M-MRS-1-2-3-PRM-0131-V1.0/DOCUMENT/ESA_DOC/IFMS_OCCFTP_10_3_1.PDF>`_ and `mission-specific documentation <https://pdssbn.astro.umd.edu/holdings/ro-x-rsi-1_2_3-cvp1-0008-v1.0/document/rsi_doc/ros_rsi_igm_ds_3118.pdf>`_. IFMS files are text files containing the radio tracking data and station transmission frequencies. The :func:`~tudatpy.estimation.observations_setup.observations_from_ifms_files` function is used to load an IFMS file, store the relevant data in an :class:`~tudatpy.estimation.observations.ObservationCollection`, and update the transmitting station's transmission frequencies.
* Fdets files: Data files containing open-loop Doppler data produced using the Planetary Radio Interferometry and Doppler Experiment (PRIDE), see :func:`data processing description <https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/high-spectral-resolution-multitone-spacecraft-doppler-tracking-software-algorithms-and-implementations/629B2AC90B0C233C0B43D0AEF897FCA6>`_. These data files contain frequency observations, typically from a spacecraft downlink that is retransmitted from an uplink. The Fdets files only contain the frequency detections. To use these data in an orbit estimation, the transmitter frequency needs to be defined, or loaded from another data source (such as ODF or IFMS files).

We also provide a generic interface for loading simple text files with data (see :func:`~tudatpy.estimation.observations_setup.observations_wrapper.create_tracking_txtfile_observation_collection`) and for creating 'pseudo-observations' (Cartesian position observations from existing ephemerides, see :func:`~tudatpy.estimation.observations_setup.observations_wrapper.create_pseudo_observations_and_models`). Finally, the :func:`~tudatpy.estimation.observations_setup.observations_wrapper.simulate_observations` function can be used to simulate observations inside Tudat.

In each of the above cases, the functionality in this module creates an :class:`~tudatpy.estimation.observations.ObservationCollection` object, which can then be used as input to perform an estimation (see :ref:`estimation_analysis` module).

Functions
---------
.. currentmodule:: tudatpy.estimation.observations_setup.observations_wrapper

.. autosummary::

   process_odf_data_single_file

   process_odf_data_multiple_files

   set_odf_information_in_bodies

   create_odf_observed_observation_collection

   observations_from_odf_files

   observations_from_ifms_files

   observations_from_multi_station_ifms_files

   observations_from_fdets_files

   create_compressed_doppler_collection

   create_tracking_txtfile_observation_collection

   create_pseudo_observations_and_models

   set_existing_observations

   simulate_observations

   single_type_observation_collection


.. autofunction:: tudatpy.estimation.observations_setup.observations_wrapper.process_odf_data_single_file
   
.. autofunction:: tudatpy.estimation.observations_setup.observations_wrapper.process_odf_data_multiple_files

.. autofunction:: tudatpy.estimation.observations_setup.observations_wrapper.set_odf_information_in_bodies

.. autofunction:: tudatpy.estimation.observations_setup.observations_wrapper.create_odf_observed_observation_collection

.. autofunction:: tudatpy.estimation.observations_setup.observations_wrapper.observations_from_odf_files

.. autofunction:: tudatpy.estimation.observations_setup.observations_wrapper.observations_from_ifms_files

.. autofunction:: tudatpy.estimation.observations_setup.observations_wrapper.observations_from_multi_station_ifms_files

.. autofunction:: tudatpy.estimation.observations_setup.observations_wrapper.observations_from_fdets_files

.. autofunction:: tudatpy.estimation.observations_setup.observations_wrapper.create_compressed_doppler_collection

.. autofunction:: tudatpy.estimation.observations_setup.observations_wrapper.create_tracking_txtfile_observation_collection

.. autofunction:: tudatpy.estimation.observations_setup.observations_wrapper.create_pseudo_observations_and_models

.. autofunction:: tudatpy.estimation.observations_setup.observations_wrapper.set_existing_observations

.. autofunction:: tudatpy.estimation.observations_setup.observations_wrapper.simulate_observations

.. autofunction:: tudatpy.estimation.observations_setup.observations_wrapper.single_type_observation_collection

Classes
-------
.. currentmodule:: tudatpy.estimation.observations_setup.observations_wrapper

.. autosummary::

   ProcessedOdfFileContents

.. autoclass:: tudatpy.estimation.observations_setup.observations_wrapper.ProcessedOdfFileContents
   :members:
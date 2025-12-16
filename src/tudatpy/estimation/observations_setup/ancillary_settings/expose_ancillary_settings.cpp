/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_ancillary_settings.h"
#include <pybind11/functional.h>
#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"

namespace tom = tudat::observation_models;
namespace tss = tudat::simulation_setup;

namespace tudat
{

namespace simulation_setup
{

void addAncillarySettingsToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >& observationSimulationSettings,
        const std::shared_ptr< tom::ObservationAncillarySimulationSettings >& ancillarySettings,
        const tom::ObservableType observableType )
{
    tss::addAncillarySettingsToObservationSimulationSettings< TIME_TYPE, const tom::ObservableType >(
            observationSimulationSettings, ancillarySettings, observableType );
}

void addAncillarySettingsToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >& observationSimulationSettings,
        const std::shared_ptr< tom::ObservationAncillarySimulationSettings >& ancillarySettings,
        const tom::ObservableType observableType,
        const tom::LinkDefinition& linkEnds )
{
    tss::addAncillarySettingsToObservationSimulationSettings< TIME_TYPE, const tom::ObservableType, const tom::LinkDefinition& >(
            observationSimulationSettings, ancillarySettings, observableType, linkEnds );
}

}  // namespace simulation_setup

}  // namespace tudat

namespace tudatpy
{
namespace estimation
{
namespace observations_setup
{

namespace ancillary_settings
{

void expose_ancillary_settings( py::module& m )
{
    py::enum_< tom::ObservationAncillarySimulationVariable >( m,
                                                               "ObservationAncillarySimulationVariable",
                                                               R"doc(

      Enumeration of observation ancillary variable types.

      )doc" )
            .value( "link_ends_delays",
                    tom::ObservationAncillarySimulationVariable::link_ends_delays,
                    R"doc(
                    Retransmission delays at the retransmitter link ends (in seconds), typically for an n-way range or Doppler observation.
                    For a set of link ends consisting of :math:`N` one-way link ends (for instance ``transmitter``->``retransmitter``->``receiver``
                    for :math:`N=2`, this ancillary setting is a list of length :math:`N-1` representing the time in seconds between
                    signal reception and subsequent retransmission/reflection at a link end. This ancillary setting is retrieved and set using the
                    :attr:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings.get_float_list_settings` and
                    :attr:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings.set_float_list_settings`
                    )doc" )
            .value( "doppler_integration_time",
                    tom::ObservationAncillarySimulationVariable::doppler_integration_time,
                    R"doc(
                    Time interval :math:`\Delta t` in seconds over which averaged Doppler observables (see :ref:`model_settings` for ``_averaged`` observation models)
                    are averaged. This quantity is also termed the count time. This ancillary setting is retrieved and set using the
                    :attr:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings.get_float_settings` and
                    :attr:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings.set_float_settings`
                    )doc" )
            .value( "doppler_reference_frequency",
                    tom::ObservationAncillarySimulationVariable::doppler_reference_frequency,
                    R"doc(
                    Reference frequency :math:`f_{\text{ref}}` w.r.t. which the Doppler observable is computed for the
                    :func:`~tudatpy.estimation.observable_models_setup.model_settings.dsn_n_way_doppler_averaged` observation model.
                    This ancillary setting is retrieved and set using the
                    :attr:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings.get_float_settings` and
                    :attr:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings.set_float_settings`
                    )doc" )
            .value( "frequency_bands",
                    tom::ObservationAncillarySimulationVariable::frequency_bands,
                    R"doc(
                    Frequency bands for the up and down-link of a radio observable, used to compute the turnaround ratio :math:`M_{2}` at the
                    retransmitter (so that the received frequency :math:`f_{1}` and the retransmitted frequency :math:`f_{2}` are related as
                    :math:`f_{2}=M_{2}\cdot f_{1}`.
                    This ancillary setting is retrieved and set using the
                    :attr:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings.get_float_list_settings` and
                    :attr:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings.set_float_list_settings`
                    the mapping of the :class:`~tudatpy.estimation.observations_setup.ancillary_settings.FrequencyBands`
                    enum to integers (with s-band equal to 1, x-band to 1, ku band to 2, ka-band to 3). So, for an s-band up and x-band downlink, this
                    ancillary setting gets the value ``[0, 1]``
                    )doc" )
            .value( "reception_reference_frequency_band",
                    tom::ObservationAncillarySimulationVariable::reception_reference_frequency_band,
                    R"doc(
                    Receiver reference frequency  band w.r.t. which the reference turnaround ratio :math:`M_{2,R}` is computed for the
                    :func:`~tudatpy.estimation.observable_models_setup.model_settings.dsn_n_way_doppler_averaged` observation model.
                    This ancillary setting is retrieved and set using the
                    :attr:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings.get_float_settings` and
                    :attr:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings.set_float_settings`
                    )doc" )
            .value( "sequential_range_lowest_ranging_component",
                    tom::ObservationAncillarySimulationVariable::sequential_range_lowest_ranging_component,
                    R"doc(
                    Lowest sequential ranging component :math:`n` used for the
                    :func:`~tudatpy.estimation.observable_models_setup.model_settings.dsn_n_way_range` observation model.
                    This ancillary setting is retrieved and set using the
                    :attr:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings.get_float_settings` and
                    :attr:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings.set_float_settings`
                    )doc" )
            .export_values( );

    py::class_< tom::ObservationAncillarySimulationSettings, std::shared_ptr< tom::ObservationAncillarySimulationSettings > >(
            m,
            "ObservationAncillarySimulationSettings",
            R"doc(

    Class for holding ancillary settings for observation simulation (see module level documentation for typical usage and creation).

    This class holds both single-valued (float) and multi-valued (list of floats) ancillary settings

      )doc" )
            .def( py::init<>( ),
                  R"doc(

                 Create an empty ancillary settings object

                 Returns
                 -------
                 estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings

                 )doc" )
            .def( "set_float_settings",
                  &tudat::observation_models::ObservationAncillarySimulationSettings::setAncillaryDoubleData,
                  py::arg( "variable" ),
                  py::arg( "value" ),
                  R"doc(

                Function to set a multi-valued ancillary setting in this object

                Parameters
                ----------
                setting_type : ObservationAncillarySimulationVariable
                   Type of the setting for which the value is to be set

                value : list[float]
                   List of value that define the provided setting type

                )doc" )
            .def( "set_float_list_settings",
                  &tudat::observation_models::ObservationAncillarySimulationSettings::setAncillaryDoubleVectorData,
                  py::arg( "variable" ),
                  py::arg( "value" ),
                  R"doc(

                Function to set a single-valued ancillary setting value in this object

                Parameters
                ----------
                setting_type : ObservationAncillarySimulationVariable
                   Type of the setting for which the value is to be set

                value : float
                   Value for the setting

                )doc" )
            .def( "set_intermediate_double_data",
                  &tudat::observation_models::ObservationAncillarySimulationSettings::setIntermediateDoubleData,
                  py::arg( "variable" ),
                  py::arg( "value" ) )
            .def( "get_float_settings",
                  &tom::ObservationAncillarySimulationSettings::getAncillaryDoubleData,
                  py::arg( "setting_type" ),
                  py::arg( "throw_exception" ) = true,
                  R"doc(

         Function to retrieve a single-valued (``float``) ancillary setting value

         Parameters
         ----------
         setting_type : ObservationAncillarySimulationVariable
             Type of the setting for which the value is to be returned

         throw_exception : bool, default = false
             Boolean defining whether to throw an exception if the requested setting does not exist, or does not exist as a floating point value.

         Returns
         -------
         float
             Value of the requested ancillary variable (or NaN if it does not exist and ``throw_exception`` is ``false``)

     )doc" )
            .def( "get_float_list_settings",
                  &tom::ObservationAncillarySimulationSettings::getAncillaryDoubleVectorData,
                  py::arg( "setting_type" ),
                  py::arg( "throw_exception" ) = true,
                  R"doc(

         Function to retrieve a multi-valued (``list`` of ``float``) ancillary setting value

         Parameters
         ----------
         setting_type : ObservationAncillarySimulationVariable
             Type of the setting for which the value is to be returned

         throw_exception : bool, default = false
             Boolean defining whether to throw an exception if the requested setting does not exist, or does not exist as list of floating point values.

         Returns
         -------
         list[ float ]
             Value of the requested ancillary variable (or empty list if it does not exist and ``throw_exception`` is ``false``)

         Examples
         --------
         .. code-block:: python

             # Code snippet to show how to retrieve ObservationAncillarySimulationSettings info
             # using the ObservationAncillarySimulationSettings.get_float_settings function

             from tudatpy.estimation.observations_setup import ancillary_settings

             # Create Ancillary Settings
             n_way_range_ancillary_settings = ancillary_settings.n_way_range_ancillary_settings(frequency_bands=[ancillary_settings.FrequencyBands.x_band])

             # Verify that we indeed added Frequency Bands as Ancillary Simulation Variables, using n_way_range_ancillary_settings.get_float_list_settings
             list_num = n_way_range_ancillary_settings.get_float_list_settings(ancillary_settings.ObservationAncillarySimulationVariable.frequency_bands)
             for num in list_num:
                 name = ancillary_settings.ObservationAncillarySimulationVariable(int(num)).name
                 print(f'Ancillary Simulation Variable(s): {name}, corresponding to enumeration object n. {int(num)} of the ObservationAncillarySimulationVariable Enumeration')




     )doc" );

    py::enum_< tudat::observation_models::ObservationIntermediateSimulationVariable >( m,
                                                                                       "ObservationIntermediateSimulationVariable",
                                                                                       R"doc(
        Enumeration of observation intermediate variable types.

        This enum lists variables that are computed during the observation simulation process and can be stored for later analysis.
        )doc" )
            .value( "transmitter_frequency_intermediate",
                    tudat::observation_models::ObservationIntermediateSimulationVariable::transmitter_frequency_intermediate )
            .value( "received_frequency_intermediate",
                    tudat::observation_models::ObservationIntermediateSimulationVariable::received_frequency_intermediate )
            .export_values( );

    m.def( "doppler_ancillary_settings",
           &tom::getAveragedDopplerAncillarySettings,
           py::arg( "integration_time" ) = 60.0,
           R"doc(

 Function for creating ancillary settings for averaged Doppler observable.

 Function for creating ancillary settings for an averaged Doppler observable. Specifically, this
 function can be used to create settings for the integration time :math:`\Delta t` of the observable
 (see :attr:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationVariable.doppler_integration_time`).
 Note: in case no retransmission
 delays (or other additional ancillary settings) are to be defined, this setting may be used for one-, two-, or N-way
 averaged Doppler.

 Parameters
 ----------
 integration_time : float, default = 60.0
     Integration time that is to be used for the averaged Doppler observable
 Returns
 -------
 :class:`ObservationAncillarySimulationSettings`
     Instance of the :class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings` with the required settings.






     )doc" );

    m.def( "two_way_range_ancillary_settings",
           &tom::getTwoWayRangeAncillarySettings,
           py::arg( "retransmission_delay" ) = 0.0,
           // py::arg("frequency_band") =
           // tom::FrequencyBands::x_band,
           R"doc(

 Function for creating ancillary settings for two-way range observable.

 Function for creating ancillary settings for a two-way range observable. Specifically, this
 function can be used to create settings for the retransmission delay of the observable (see :attr:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationVariable.link_ends_delays`).
 Note: this function is provided for convenience, and is equivalent to calling :func:`~tudatpy.estimation.observations_setup.ancillary_settings.n_way_range_ancillary_settings`
 with a single retransmission delay.


 Parameters
 ----------
 retransmission_delay : float, default = 0.0
     Retransmission delay that is to be applied to the simulation of the two-way observable
 Returns
 -------
 :class:`ObservationAncillarySimulationSettings`
     Instance of the :class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings` with the required settings.






     )doc" );

    m.def( "two_way_doppler_ancillary_settings",
           &tom::getTwoWayAveragedDopplerAncillarySettings,
           py::arg( "integration_time" ) = 60.0,
           py::arg( "retransmission_delay" ) = 0.0,
           R"doc(

 Function for creating ancillary settings for two-way averaged Doppler observable.

 Function for creating ancillary settings for a two-way range observable. Specifically, this
 function can be used to create settings for the retransmission delay of the observable  (see :attr:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationVariable.link_ends_delays`).
 Note:
 this function is provided for convenience, and is equivalent to calling :func:`~tudatpy.estimation.observations_setup.ancillary_settings.n_way_doppler_ancillary_settings`
 with a single retransmission delay.


 Parameters
 ----------
 integration_time : float, default = 60.0
     Integration time that is to be used for the averaged Doppler observable
 retransmission_delay : float, default = 0.0
     Retransmission delay that is to be applied to the simulation of the two-way observable
 Returns
 -------
 :class:`ObservationAncillarySimulationSettings`
     Instance of the :class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings` with the required settings.






     )doc" );

    m.def( "n_way_range_ancillary_settings",
           &tom::getNWayRangeAncillarySettings,
           py::arg( "link_end_delays" ) = std::vector< double >( ),
           py::arg( "frequency_bands" ) = std::vector< tom::FrequencyBands >( ),
           R"doc(

 Function for creating ancillary settings for n-way range observable.

 Function for creating ancillary settings for a n-way range observable. Specifically, this
 function can be used to create settings for the retransmission delays of the observable (see :attr:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationVariable.link_ends_delays`), for each of the retransmitters.


 Parameters
 ----------
 retransmission_delays : list[ float ], default = None
     Retransmission delays that are to be applied to the simulation of the n-way observable. If kept empty, this results in 0 retransmission delay at each retransmitter. If defined, this list must be the same length as the number of retransmitters, and the :math:`i^{th}` entry contains the retransmission delay of the :math:`i^{th}` retrasmitter
 Returns
 -------
 :class:`ObservationAncillarySimulationSettings`
     Instance of the :class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings` with the required settings.






     )doc" );

    m.def( "n_way_doppler_ancillary_settings",
           &tom::getNWayAveragedDopplerAncillarySettings,
           py::arg( "integration_time" ) = 60.0,
           py::arg( "link_end_delays" ) = std::vector< double >( ),
           py::arg( "frequency_bands" ) = std::vector< tom::FrequencyBands >( ),
           R"doc(

 Function for creating ancillary settings for n-way averaged Doppler observable.

 Function for creating ancillary settings for a n-way averaged Doppler observable. Specifically, this
 function can be used to create settings for the integration time of the observable (see :attr:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationVariable.doppler_integration_time`),
 and the  retransmission delays for each of the retransmitters (see :attr:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationVariable.link_ends_delays`).


 Parameters
 ----------
 integration_time : float, default = 60.0
     Integration time that is to be used for the averaged Doppler observable
 retransmission_delays : list[ float ], default = None
     Retransmission delays that are to be applied to the simulation of the n-way observable. If kept empty, this results in 0 retransmission delay at each retransmitter. If defined, this list must be the same length as the number of retransmitters, and the :math:`i^{th}` entry contains the retransmission delay of the :math:`i^{th}` retrasmitter
 Returns
 -------
 :class:`ObservationAncillarySimulationSettings`
     Instance of the :class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings` with the required settings.






     )doc" );

    m.def( "dsn_n_way_doppler_ancillary_settings",
           &tom::getDsnNWayAveragedDopplerAncillarySettings,
           py::arg( "frequency_bands" ),
           py::arg( "reference_frequency_band" ),
           py::arg( "reference_frequency" ),
           py::arg( "integration_time" ) = 60.0,
           py::arg( "link_end_delays" ) = std::vector< double >( ),
           R"doc(

        Function for creating ancillary settings for DSN n-way averaged Doppler observable (:func:`~tudatpy.estimation.observable_models_setup.model_settings.dsn_n_way_doppler_averaged`).

        Parameters
        ----------
        frequency_bands : list[tudatpy.estimation.observations_setup.ancillary_settings.FrequencyBands]
            List of frequency bands for each link.
        reference_frequency_band : tudatpy.estimation.observations_setup.ancillary_settings.FrequencyBands
            Reference frequency band at reception.
        reference_frequency : float
            Reference frequency for Doppler calculation.
        integration_time : float, default = 60.0
            Integration time for the averaged Doppler observable.
        link_end_delays : list[float], default = []
            Retransmission delays at each retransmitter.

        Returns
        -------
        tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings
            Instance of the ObservationAncillarySimulationSettings with the required settings.
        )doc" );

    m.def( "dsn_n_way_range_ancillary_settings",
           &tom::getDsnNWayRangeAncillarySettings,
           py::arg( "frequency_bands" ),
           py::arg( "lowest_ranging_component" ),
           py::arg( "link_end_delays" ) = std::vector< double >( ),
           R"doc(
        Function for creating ancillary settings for DSN n-way range observable (:func:`~tudatpy.estimation.observable_models_setup.model_settings.dsn_n_way_range`).

        Parameters
        ----------
        frequency_bands : list[tudatpy.estimation.observations_setup.ancillary_settings.FrequencyBands]
            List of frequency bands for each link.
        lowest_ranging_component : float
            Lowest ranging component for sequential range.
        link_end_delays : list[float], default = []
            Retransmission delays at each retransmitter.

        Returns
        -------
        tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings
            Instance of the ObservationAncillarySimulationSettings with the required settings.
        )doc" );

    m.def( "doppler_measured_frequency_ancillary_settings",
           &tom::getDopplerMeasuredFrequencyAncillarySettings,
           py::arg( "frequency_bands" ),
           R"doc(

        Function for creating ancillary settings for two-way Doppler frequency observable (:func:`~tudatpy.estimation.observable_models_setup.model_settings.two_way_doppler_instantaneous_frequency`).

        Parameters
        ----------
        frequency_bands : list[tudatpy.estimation.observations_setup.ancillary_settings.FrequencyBands]
            List of frequency bands for each link.

        Returns
        -------
        tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings
            Instance of the ObservationAncillarySimulationSettings with the required settings.
        )doc" );

    m.def( "add_ancillary_settings_to_observable",
           py::overload_cast< const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const std::shared_ptr< tom::ObservationAncillarySimulationSettings >&,
                              const tom::ObservableType >( &tss::addAncillarySettingsToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "ancillary_settings" ),
           py::arg( "observable_type" ),
           R"doc(

        Add ancillary settings to all observation simulation settings of a given observable type.

        Parameters
        ----------
        observation_simulation_settings_list : list[tudatpy.estimation.observations_setup.ObservationSimulationSettings]
            List of observation simulation settings to modify.
        ancillary_settings : tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings
            Ancillary settings to add.
        observable_type : tudatpy.kernel.astro.ObservableType
            Observable type for which to add the ancillary settings.
        )doc" );

    m.def( "add_ancillary_settings_to_observable_for_link_ends",
           py::overload_cast< const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const std::shared_ptr< tom::ObservationAncillarySimulationSettings >&,
                              const tom::ObservableType,
                              const tom::LinkDefinition& >( &tss::addAncillarySettingsToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "ancillary_settings" ),
           py::arg( "observable_type" ),
           py::arg( "link_ends" ),
           R"doc(
        Add ancillary settings to observation simulation settings for a specific observable type and link ends.

        Parameters
        ----------
        observation_simulation_settings_list : list[tudatpy.estimation.observations_setup.ObservationSimulationSettings]
            List of observation simulation settings to modify.
        ancillary_settings : tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings
            Ancillary settings to add.
        observable_type : tudatpy.kernel.astro.ObservableType
            Observable type for which to add the ancillary settings.
        link_ends : tudatpy.kernel.astro.LinkDefinition
            Link ends for which to add the ancillary settings.
        )doc" );

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // FREQUENCIES
    /////////////////////////////////////////////////////////////////////////////////////////////////

    py::enum_< tom::FrequencyBands >( m, "FrequencyBands", R"doc(
        Enumeration of frequency bands.

        This enum lists common frequency bands used in deep space navigation.
        )doc" )
            .value( "s_band", tom::FrequencyBands::s_band )
            .value( "x_band", tom::FrequencyBands::x_band )
            .value( "ka_band", tom::FrequencyBands::ka_band )
            .value( "ku_band", tom::FrequencyBands::ku_band );

    m.def( "dsn_default_turnaround_ratios",
           &tom::getDsnDefaultTurnaroundRatios,
           py::arg( "uplink_band" ),
           py::arg( "downlink_band" ),
           R"doc(

        Get the default DSN turnaround ratio for given uplink and downlink bands.

        Parameters
        ----------
        uplink_band : tudatpy.estimation.observations_setup.ancillary_settings.FrequencyBands
            Uplink frequency band.
        downlink_band : tudatpy.estimation.observations_setup.ancillary_settings.FrequencyBands
            Downlink frequency band.

        Returns
        -------
        float
            The turnaround ratio.
        )doc" );

    m.def( "cassini_turnaround_ratios",
           &tom::getCassiniTurnaroundRatio,
           py::arg( "uplink_band" ),
           py::arg( "downlink_band" ),
           R"doc(
        Get the Cassini turnaround ratio for given uplink and downlink bands.

        This function returns the specific Ka-band turnaround ratio for Cassini if both bands are Ka-band,
        otherwise it falls back to the DSN default turnaround ratios.

        Parameters
        ----------
        uplink_band : tudatpy.estimation.observations_setup.ancillary_settings.FrequencyBands
            Uplink frequency band.
        downlink_band : tudatpy.estimation.observations_setup.ancillary_settings.FrequencyBands
            Downlink frequency band.

        Returns
        -------
        float
            The turnaround ratio.
        )doc" );
}

}  // namespace ancillary_settings
}  // namespace observations_setup
}  // namespace estimation
}  // namespace tudatpy
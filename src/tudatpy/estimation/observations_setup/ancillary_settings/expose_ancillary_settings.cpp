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

void addAncilliarySettingsToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >& observationSimulationSettings,
        const std::shared_ptr< tom::ObservationAncilliarySimulationSettings >& ancilliarySettings,
        const tom::ObservableType observableType )
{
    tss::addAncilliarySettingsToObservationSimulationSettings< TIME_TYPE, const tom::ObservableType >(
            observationSimulationSettings, ancilliarySettings, observableType );
}

void addAncilliarySettingsToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >& observationSimulationSettings,
        const std::shared_ptr< tom::ObservationAncilliarySimulationSettings >& ancilliarySettings,
        const tom::ObservableType observableType,
        const tom::LinkDefinition& linkEnds )
{
    tss::addAncilliarySettingsToObservationSimulationSettings< TIME_TYPE, const tom::ObservableType, const tom::LinkDefinition& >(
            observationSimulationSettings, ancilliarySettings, observableType, linkEnds );
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
    py::enum_< tom::ObservationAncilliarySimulationVariable >( m,
                                                               "ObservationAncilliarySimulationVariable",
                                                               R"doc(

Enumeration of observation ancillary variable types.

Examples
--------

.. code-block:: python

    # Code snippet to print all available Observation Ancillary Variable Types
    from tudatpy.estimation.observations_setup import ancillary_settings

    num_observation_ancillary_variable_types = len(ancillary_settings.ObservationAncilliarySimulationVariable.__members__)
    print(f'The length of all available Tudatpy  Observation Ancillary Variable Types is: {num_observation_ancillary_variable_types}')

    # Print all Observation Ancillary Variable Types using the "name" property
    for i in range(num_observation_ancillary_variable_types):
        print(i, ancillary_settings.ObservationAncilliarySimulationVariable(i).name)



      )doc" )
            .value( "link_ends_delays", tom::ObservationAncilliarySimulationVariable::link_ends_delays )
            .value( "doppler_integration_time", tom::ObservationAncilliarySimulationVariable::doppler_integration_time )
            .value( "doppler_reference_frequency", tom::ObservationAncilliarySimulationVariable::doppler_reference_frequency )
            .value( "frequency_bands", tom::ObservationAncilliarySimulationVariable::frequency_bands )
            .value( "reception_reference_frequency_band", tom::ObservationAncilliarySimulationVariable::reception_reference_frequency_band )
            .value( "sequential_range_lowest_ranging_component",
                    tom::ObservationAncilliarySimulationVariable::sequential_range_lowest_ranging_component )
            .value( "range_conversion_factor", tom::ObservationAncilliarySimulationVariable::range_conversion_factor )
            .export_values( );

    py::class_< tom::ObservationAncilliarySimulationSettings, std::shared_ptr< tom::ObservationAncilliarySimulationSettings > >(
            m,
            "ObservationAncilliarySimulationSettings",
            R"doc(

         Base class for holding ancilliary settings for observation simulation.

         Base class for holding ancilliary settings for observation simulation.
         The user can create instances of this class via the :func:`~estimation.observations_setup.observations_dependent_variables.elevation_angle_dependent_variable` function.

         Examples
         --------
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

             # [OPTIONAL] Verify that we indeed added Frequency Bands as Ancillary Simulation Variables for the n_way_range_ancillary_settings.
             list_num = n_way_range_ancillary_settings.get_float_list_settings(ancillary_settings.ObservationAncilliarySimulationVariable.frequency_bands)
             for num in list_num:
                 name = ancillary_settings.ObservationAncilliarySimulationVariable(int(num)).name
                 print(f'Ancillary Simulation Variable(s): {name}, corresponding to enumeration object n. {int(num)} of the ObservationAncilliarySimulationVariable Enumeration')



      )doc" )
            .def( "set_intermediate_double_data",
                  &tudat::observation_models::ObservationAncilliarySimulationSettings::setIntermediateDoubleData,
                  py::arg( "variable" ),
                  py::arg( "value" ) )
            .def( "get_float_settings",
                  &tom::ObservationAncilliarySimulationSettings::getAncilliaryDoubleData,
                  py::arg( "setting_type" ),
                  py::arg( "throw_exception" ) = true,
                  R"doc(


         Parameters
         ----------
         setting_type : ObservationAncilliarySimulationVariable
             Type of the setting for which the value is to be returned

         throw_exception : bool, default = false
             Boolean defining whether to throw an exception if the requested setting does not exist, or does not exist as a floating point value.

         Returns
         -------
         float
             Value of the requested ancilliary variable (or NaN if it does not exist)

     )doc" )
            .def( "get_float_list_settings",
                  &tom::ObservationAncilliarySimulationSettings::getAncilliaryDoubleVectorData,
                  py::arg( "setting_type" ),
                  py::arg( "throw_exception" ) = true,
                  R"doc(


         Parameters
         ----------
         setting_type : ObservationAncilliarySimulationVariable
             Type of the setting for which the value is to be returned

         throw_exception : bool, default = false
             Boolean defining whether to throw an exception if the requested setting does not exist, or does not exist as list of floating point values.

         Returns
         -------
         list[ float ]
             Value of the requested ancilliary variable (or empty list if it does not exist)

         Examples
         --------
         .. code-block:: python

             # Code snippet to show how to retrieve ObservationAncillarySimulationSettings info
             # using the ObservationAncilliarySimulationSettings.get_float_settings function

             from tudatpy.estimation.observations_setup import ancillary_settings

             # Create Ancillary Settings
             n_way_range_ancillary_settings = ancillary_settings.n_way_range_ancilliary_settings(frequency_bands=[ancillary_settings.FrequencyBands.x_band])

             # Verify that we indeed added Frequency Bands as Ancillary Simulation Variables, using n_way_range_ancillary_settings.get_float_list_settings
             list_num = n_way_range_ancillary_settings.get_float_list_settings(ancillary_settings.ObservationAncilliarySimulationVariable.frequency_bands)
             for num in list_num:
                 name = ancillary_settings.ObservationAncilliarySimulationVariable(int(num)).name
                 print(f'Ancillary Simulation Variable(s): {name}, corresponding to enumeration object n. {int(num)} of the ObservationAncilliarySimulationVariable Enumeration')




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

    m.def( "doppler_ancilliary_settings",
           &tom::getAveragedDopplerAncilliarySettings,
           py::arg( "integration_time" ) = 60.0,
           R"doc(

 Function for creating ancilliary settings for averaged Doppler observable.

 Function for creating ancilliary settings for an averaged Doppler observable. Specifically, this
 function can be used to create settings for the integration time of the observable. Note: in case no retransmission
 delays (or other additional ancilliary settings) are to be defined, this setting may be used for one-, two-, or N-way
 averaged Doppler.


 Parameters
 ----------
 integration_time : float, default = 60.0
     Integration time that is to be used for the averaged Doppler observable
 Returns
 -------
 :class:`ObservationAncilliarySimulationSettings`
     Instance of the :class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings` with the required settings.






     )doc" );

    m.def( "two_way_range_ancilliary_settings",
           &tom::getTwoWayRangeAncilliarySettings,
           py::arg( "retransmission_delay" ) = 0.0,
           // py::arg("frequency_band") =
           // tom::FrequencyBands::x_band,
           R"doc(

 Function for creating ancilliary settings for two-way range observable.

 Function for creating ancilliary settings for a two-way range observable. Specifically, this
 function can be used to create settings for the retransmission delay of the observable. NOTE:
 this function is provided for convenience, and is equivalent to calling :func:`~tudatpy.estimation.observations_setup.ancillary_settings.n_way_range_ancilliary_settings`
 with a single retransmission delay.


 Parameters
 ----------
 retransmission_delay : float, default = 0.0
     Retransmission delay that is to be applied to the simulation of the two-way observable
 Returns
 -------
 :class:`ObservationAncilliarySimulationSettings`
     Instance of the :class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings` with the required settings.






     )doc" );

    m.def( "two_way_doppler_ancilliary_settings",
           &tom::getTwoWayAveragedDopplerAncilliarySettings,
           py::arg( "integration_time" ) = 60.0,
           py::arg( "retransmission_delay" ) = 0.0,
           R"doc(

 Function for creating ancilliary settings for two-way averaged Doppler observable.

 Function for creating ancilliary settings for a two-way range observable. Specifically, this
 function can be used to create settings for the retransmission delay of the observable.  NOTE:
 this function is provided for convenience, and is equivalent to calling :func:`~tudatpy.estimation.observations_setup.ancillary_settings.n_way_doppler_ancilliary_settings`
 with a single retransmission delay.


 Parameters
 ----------
 integration_time : float, default = 60.0
     Integration time that is to be used for the averaged Doppler observable
 retransmission_delay : float, default = 0.0
     Retransmission delay that is to be applied to the simulation of the two-way observable
 Returns
 -------
 :class:`ObservationAncilliarySimulationSettings`
     Instance of the :class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings` with the required settings.






     )doc" );

    m.def( "n_way_range_ancilliary_settings",
           &tom::getNWayRangeAncilliarySettings,
           py::arg( "link_end_delays" ) = std::vector< double >( ),
           py::arg( "frequency_bands" ) = std::vector< tom::FrequencyBands >( ),
           R"doc(

 Function for creating ancilliary settings for n-way range observable.

 Function for creating ancilliary settings for a n-way range observable. Specifically, this
 function can be used to create settings for the retransmission delays of the observable, for each of the retransmitters.


 Parameters
 ----------
 retransmission_delays : list[ float ], default = None
     Retransmission delays that are to be applied to the simulation of the n-way observable. If kept empty, this results in 0 retransmission delay at each retransmitter. If defined, this list must be the same length as the number of retransmitters, and the :math:`i^{th}` entry contains the retransmission delay of the :math:`i^{th}` retrasmitter
 Returns
 -------
 :class:`ObservationAncilliarySimulationSettings`
     Instance of the :class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings` with the required settings.






     )doc" );

    m.def( "n_way_doppler_ancilliary_settings",
           &tom::getNWayAveragedDopplerAncilliarySettings,
           py::arg( "integration_time" ) = 60.0,
           py::arg( "link_end_delays" ) = std::vector< double >( ),
           py::arg( "frequency_bands" ) = std::vector< tom::FrequencyBands >( ),
           R"doc(

 Function for creating ancilliary settings for n-way averaged Doppler observable.

 Function for creating ancilliary settings for a n-way averaged Doppler observable. Specifically, this
 function can be used to create settings for the integration time of the observable, and the  retransmission delays for each of the retransmitters.


 Parameters
 ----------
 integration_time : float, default = 60.0
     Integration time that is to be used for the averaged Doppler observable
 retransmission_delays : list[ float ], default = None
     Retransmission delays that are to be applied to the simulation of the n-way observable. If kept empty, this results in 0 retransmission delay at each retransmitter. If defined, this list must be the same length as the number of retransmitters, and the :math:`i^{th}` entry contains the retransmission delay of the :math:`i^{th}` retrasmitter
 Returns
 -------
 :class:`ObservationAncilliarySimulationSettings`
     Instance of the :class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings` with the required settings.






     )doc" );

    m.def( "dsn_n_way_doppler_ancilliary_settings",
           &tom::getDsnNWayAveragedDopplerAncillarySettings,
           py::arg( "frequency_bands" ),
           py::arg( "reference_frequency_band" ),
           py::arg( "reference_frequency" ),
           py::arg( "integration_time" ) = 60.0,
           py::arg( "link_end_delays" ) = std::vector< double >( ),
           R"doc(
        Function for creating ancilliary settings for DSN n-way averaged Doppler observable.

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
        tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings
            Instance of the ObservationAncilliarySimulationSettings with the required settings.
        )doc" );

    m.def( "dsn_n_way_range_ancilliary_settings",
           &tom::getDsnNWayRangeAncillarySettings,
           py::arg( "frequency_bands" ),
           py::arg( "lowest_ranging_component" ),
           py::arg( "link_end_delays" ) = std::vector< double >( ),
           R"doc(
        Function for creating ancilliary settings for DSN n-way range observable.

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
        tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings
            Instance of the ObservationAncilliarySimulationSettings with the required settings.
        )doc" );

    m.def( "doppler_measured_frequency_ancillary_settings",
           &tom::getDopplerMeasuredFrequencyAncilliarySettings,
           py::arg( "frequency_bands" ),
           R"doc(
        Function for creating ancilliary settings for Doppler observable with measured frequency.

        Parameters
        ----------
        frequency_bands : list[tudatpy.estimation.observations_setup.ancillary_settings.FrequencyBands]
            List of frequency bands for each link.

        Returns
        -------
        tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings
            Instance of the ObservationAncilliarySimulationSettings with the required settings.
        )doc" );

    m.def( "add_ancilliary_settings_to_observable",
           py::overload_cast< const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const std::shared_ptr< tom::ObservationAncilliarySimulationSettings >&,
                              const tom::ObservableType >( &tss::addAncilliarySettingsToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "ancilliary_settings" ),
           py::arg( "observable_type" ),
           R"doc(
        Add ancillary settings to all observation simulation settings of a given observable type.

        Parameters
        ----------
        observation_simulation_settings_list : list[tudatpy.estimation.observations_setup.ObservationSimulationSettings]
            List of observation simulation settings to modify.
        ancilliary_settings : tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings
            Ancillary settings to add.
        observable_type : tudatpy.kernel.astro.ObservableType
            Observable type for which to add the ancillary settings.
        )doc" );

    m.def( "add_ancilliary_settings_to_observable_for_link_ends",
           py::overload_cast< const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const std::shared_ptr< tom::ObservationAncilliarySimulationSettings >&,
                              const tom::ObservableType,
                              const tom::LinkDefinition& >( &tss::addAncilliarySettingsToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "ancilliary_settings" ),
           py::arg( "observable_type" ),
           py::arg( "link_ends" ),
           R"doc(
        Add ancillary settings to observation simulation settings for a specific observable type and link ends.

        Parameters
        ----------
        observation_simulation_settings_list : list[tudatpy.estimation.observations_setup.ObservationSimulationSettings]
            List of observation simulation settings to modify.
        ancilliary_settings : tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings
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
/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_observation_setup.h"

#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"

#include "tudatpy/docstrings.h"

namespace tss = tudat::simulation_setup;
namespace tom = tudat::observation_models;

namespace tudat
{

namespace simulation_setup
{

void addGaussianNoiseFunctionToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< double > > >& observationSimulationSettings,
        const double observationNoiseAmplitude,
        const tom::ObservableType observableType )
{
    tss::addGaussianNoiseFunctionToObservationSimulationSettings< double, const tom::ObservableType >(
               observationSimulationSettings, observationNoiseAmplitude, observableType );
}

void addViabilityToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< double > > >& observationSimulationSettings,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >& viabilitySettingsList )
{
    tss::addViabilityToObservationSimulationSettings< double >( observationSimulationSettings, viabilitySettingsList );
}

void addViabilityToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< double > > >& observationSimulationSettings,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >& viabilitySettingsList,
        const tom::ObservableType observableType  )
{
    tss::addViabilityToObservationSimulationSettings< double, const tom::ObservableType >(
                observationSimulationSettings, viabilitySettingsList, observableType );
}

}

}

namespace tudatpy {
namespace numerical_simulation {
namespace estimation_setup {
namespace observation {

void expose_observation_setup(py::module &m) {


    // ################      Link Definition         ################

    py::enum_< tom::LinkEndType >(m, "LinkEndType",
                                  get_docstring("LinkEndType").c_str() )
            .value("unidentified_link_end", tom::LinkEndType::unidentified_link_end )
            .value("transmitter", tom::LinkEndType::transmitter )
            .value("reflector1", tom::LinkEndType::reflector1 )
            .value("retransmitter", tom::LinkEndType::retransmitter )
            .value("reflector2", tom::LinkEndType::reflector2 )
            .value("reflector3", tom::LinkEndType::reflector3 )
            .value("reflector4", tom::LinkEndType::reflector4 )
            .value("receiver", tom::LinkEndType::receiver )
            .value("observed_body", tom::LinkEndType::observed_body )
            .export_values();

    m.def("one_way_downlink_link_ends",
          &tom::getOneWayDownlinkLinkEndsList,
          py::arg("transmitter"),
          py::arg("receivers"),
          get_docstring("one_way_downlink_link_ends").c_str() );

    m.def("one_way_uplink_link_ends",
          &tom::getOneWayUplinkLinkEndsList,
          py::arg("transmitters"),
          py::arg("receiver"),
          get_docstring("one_way_uplink_link_ends").c_str() );


    // ###########      Observation Model Settings        ################

    py::enum_< tom::ObservableType >(m, "ObservableType",
                                     get_docstring("ObservableType").c_str() )
            .value("one_way_range_type", tom::ObservableType::one_way_range )
            .value("angular_position_type", tom::ObservableType::angular_position )
            .value("position_observable_type", tom::ObservableType::position_observable )
            .value("one_way_doppler_type", tom::ObservableType::one_way_doppler )
            .value("one_way_differenced_range_type", tom::ObservableType::one_way_differenced_range )
            .value("n_way_range_type", tom::ObservableType::n_way_range )
            .value("two_way_doppler_type", tom::ObservableType::two_way_doppler )
            .value("euler_angle_313_observable_type", tom::ObservableType::euler_angle_313_observable )
            .value("velocity_observable_type", tom::ObservableType::velocity_observable )
            .export_values();


    py::class_<tom::DopplerProperTimeRateSettings,
            std::shared_ptr<tom::DopplerProperTimeRateSettings>>(
                m, "DopplerProperTimeRateSettings",
                get_docstring("DopplerProperTimeRateSettings").c_str() );

    py::class_<tom::ObservationModelSettings,
            std::shared_ptr<tom::ObservationModelSettings>>(
                m, "ObservationSettings",
                get_docstring("ObservationSettings").c_str() );

    py::class_<tom::OneWayDopplerObservationSettings,
            std::shared_ptr<tom::OneWayDopplerObservationSettings>,
            tom::ObservationModelSettings >(
            m, "OneWayDopplerObservationSettings",
            get_docstring("OneWayDopplerObservationSettings").c_str() );

    m.def("one_way_range",
          &tom::oneWayRangeSettings,
          py::arg("link_ends"),
          py::arg("light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
          py::arg("bias_settings") = nullptr,
          get_docstring("one_way_range").c_str() );

    m.def("angular_position",
          &tom::angularPositionSettings,
          py::arg("link_ends"),
          py::arg("light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
          py::arg("bias_settings") = nullptr,
          get_docstring("angular_position").c_str() );

    m.def("cartesian_position",
          &tom::positionObservableSettings,
          py::arg("link_ends"),
          py::arg("bias_settings") = nullptr,
          get_docstring("cartesian_position").c_str() );

    m.def("cartesian_velocity",
          &tom::velocityObservableSettings,
          py::arg("link_ends"),
          py::arg("bias_settings") = nullptr,
          get_docstring("cartesian_velocity").c_str() );

    m.def("313_euler_angles",
          &tom::eulerAngle313ObservableSettings,
          py::arg("link_ends"),
          py::arg("bias_settings") = nullptr,
          get_docstring("313_euler_angles").c_str() );

    m.def("one_way_open_loop_doppler",
          &tom::oneWayOpenLoopDoppler,
          py::arg("link_ends"),
          py::arg("light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
          py::arg("bias_settings") = nullptr,
          py::arg("transmitter_proper_time_rate_settings") = nullptr,
          py::arg("receiver_proper_time_rate_settings") = nullptr,
          get_docstring("313_euler_angles").c_str() );

    m.def("two_way_open_loop_doppler_from_one_way_links",
          py::overload_cast<
          const std::shared_ptr< tom::OneWayDopplerObservationSettings >,
          const std::shared_ptr< tom::OneWayDopplerObservationSettings >,
          const std::shared_ptr< tom::ObservationBiasSettings > > ( &tom::twoWayOpenLoopDoppler ),
          py::arg("uplink_doppler_settings" ),
          py::arg("downlink_doppler_settings" ),
          py::arg("bias_settings") = nullptr,
          get_docstring("two_way_open_loop_doppler_from_one_way_links").c_str() );

    m.def("two_way_open_loop_doppler",
          py::overload_cast<
          const tom::LinkEnds&,
          const std::shared_ptr< tom::LightTimeCorrectionSettings >,
          const std::shared_ptr< tom::ObservationBiasSettings > >( &tom::twoWayOpenLoopDoppler ),
          py::arg("link_ends" ),
          py::arg("light_time_correction_settings" ) = nullptr,
          py::arg("bias_settings") = nullptr,
          get_docstring("two_way_open_loop_doppler").c_str());

    m.def("one_way_closed_loop_doppler",
          py::overload_cast<
          const tom::LinkEnds&,
          const double,
          const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >,
          const std::shared_ptr< tom::ObservationBiasSettings > >( &tom::oneWayClosedLoopDoppler ),
          py::arg("link_ends" ),
          py::arg("integration_time" ),
          py::arg("light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
          py::arg("bias_settings") = nullptr,
          get_docstring("one_way_closed_loop_doppler", 0).c_str() );

    m.def("one_way_closed_loop_doppler",
          py::overload_cast<
          const tom::LinkEnds&,
          const std::function< double( const double ) >,
          const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >,
          const std::shared_ptr< tom::ObservationBiasSettings > >( &tom::oneWayClosedLoopDoppler ),
          py::arg("link_ends" ),
          py::arg("integration_time_function" ),
          py::arg("light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
          py::arg("bias_settings") = nullptr,
          get_docstring("one_way_closed_loop_doppler", 1).c_str() );


    m.def("n_way_range",
          &tom::nWayRange,
          py::arg("one_way_range_settings" ),
          py::arg("bias_settings" ) = nullptr,
          py::arg("retransmission_times_function" ) = nullptr,
          get_docstring("n_way_range").c_str() );


    py::class_<tom::LightTimeCorrectionSettings,
            std::shared_ptr<tom::LightTimeCorrectionSettings>>(
                m, "LightTimeCorrectionSettings",
                get_docstring("LightTimeCorrectionSettings").c_str() );

    m.def("first_order_relativistic_light_time_correction",
          &tom::firstOrderRelativisticLightTimeCorrectionSettings,
          py::arg("perturbing_bodies"),
          get_docstring("first_order_relativistic_light_time_correction").c_str() );

    py::class_<tom::ObservationBiasSettings,
            std::shared_ptr<tom::ObservationBiasSettings>>(
                m, "ObservationBiasSettings",
                get_docstring("ObservationBiasSettings").c_str() );

    m.def("absolute_bias",
          &tom::constantAbsoluteBias,
          py::arg("bias_value"),
          get_docstring("absolute_bias").c_str() );

    m.def("relative_bias",
          &tom::constantRelativeBias,
          py::arg("bias_value"),
          get_docstring("relative_bias").c_str() );

    m.def("arcwise_absolute_bias",
          py::overload_cast<
          const std::vector< double >&,
          const std::vector< Eigen::VectorXd >&,
          const tom::LinkEndType >( &tom::arcWiseAbsoluteBias ),
          py::arg("arc_start_times" ),
          py::arg("bias_values"),
          py::arg("reference_link_end_type" ),
          get_docstring("arcwise_absolute_bias", 0).c_str() );

    m.def("arcwise_absolute_bias",
          py::overload_cast<
          const std::map< double, Eigen::VectorXd >&,
          const tom::LinkEndType >( &tom::arcWiseAbsoluteBias ),
          py::arg("bias_values_per_start_time"),
          py::arg("reference_link_end_type" ),
          get_docstring("arcwise_absolute_bias", 1).c_str() );

    m.def("arcwise_relative_bias",
          py::overload_cast<
          const std::vector< double >&,
          const std::vector< Eigen::VectorXd >&,
          const tom::LinkEndType >( &tom::arcWiseRelativeBias ),
          py::arg("arc_start_times" ),
          py::arg("bias_values"),
          py::arg("reference_link_end_type" ),
          get_docstring("arcwise_relative_bias", 0).c_str() );

    m.def("arcwise_relative_bias",
          py::overload_cast<
          const std::map< double, Eigen::VectorXd >&,
          const tom::LinkEndType >( &tom::arcWiseRelativeBias ),
          py::arg("bias_values_per_start_time"),
          py::arg("reference_link_end_type" ),
          get_docstring("arcwise_relative_bias", 1).c_str() );

    m.def("combined_bias",
          &tom::multipleObservationBiasSettings,
          py::arg("bias_list"),
          get_docstring("combined_bias").c_str() );


    // ###########    Observation Simulation Settings     #############

    py::enum_< tom::ObservationViabilityType >(m, "ObservationViabilityType",
                                               get_docstring("ObservationViabilityType").c_str() )
            .value("minimum_elevation_angle", tom::ObservationViabilityType::minimum_elevation_angle )
            .value("body_avoidance_angle", tom::ObservationViabilityType::body_avoidance_angle )
            .value("body_occultation", tom::ObservationViabilityType::body_occultation )
            .export_values();


    py::class_<tom::ObservationViabilitySettings,
            std::shared_ptr<tom::ObservationViabilitySettings>>(
                m, "ObservationViabilitySettings",
                get_docstring("ObservationViabilityType").c_str() );

    m.def("minimum_elevation_angle_viability",
          py::overload_cast<
                  const std::pair< std::string, std::string >,
                  const double>(
                  &tom::minimumElevationAngleViabilitySettings ),
            py::arg("link_end_id" ),
            py::arg("elevation_angle" ),
            get_docstring("elevation_angle_viability").c_str() );

    m.def("maximum_elevation_angle_viability",
          py::overload_cast<
                  const std::pair< std::string, std::string >,
                  const double>(
                  &tom::maximumElevationAngleViabilitySettings ),
          py::arg("link_end_id" ),
          py::arg("elevation_angle" ),
          get_docstring("elevation_angle_viability").c_str() );

    m.def("body_avoidance_viability",
          py::overload_cast<
                  const std::pair< std::string, std::string >,
                  const std::string,
                  const double >(
                  &tom::bodyAvoidanceAngleViabilitySettings ),
          py::arg("link_end_id" ),
          py::arg("body_to_avoid" ),
          py::arg("avoidance_angle"),
          get_docstring("body_avoidance_viability").c_str() );

    m.def("body_occultation_viability",
          py::overload_cast<
                  const std::pair< std::string, std::string >,
                  const std::string >(
                  &tom::bodyOccultationViabilitySettings ),
          py::arg("link_end_id" ),
          py::arg("occulting_body" ),
          get_docstring("body_occultation_viability").c_str() );

    m.def("minimum_elevation_angle_viability_list",
          py::overload_cast<
          const std::vector< std::pair< std::string, std::string > >,
          const double >(
          &tom::minimumElevationAngleViabilitySettings ),
          py::arg("link_end_ids" ),
          py::arg("elevation_angle" ),
          get_docstring("elevation_angle_viability_list").c_str() );

    m.def("maximum_elevation_angle_viability_list",
          py::overload_cast<
                  const std::vector< std::pair< std::string, std::string > >,
                  const double >(
                  &tom::maximumElevationAngleViabilitySettings ),
          py::arg("link_end_ids" ),
          py::arg("elevation_angle" ),
          get_docstring("elevation_angle_viability_list").c_str() );

    m.def("body_avoidance_viability_list",
          py::overload_cast<
          const std::vector< std::pair< std::string, std::string > >,
          const std::string,
          const double >(
          &tom::bodyAvoidanceAngleViabilitySettings ),
          py::arg("link_end_ids" ),
          py::arg("body_to_avoid" ),
          py::arg("avoidance_angle"),
          get_docstring("body_avoidance_viability_list").c_str() );

    m.def("body_occultation_viability_list",
          py::overload_cast<
          const std::pair< std::string, std::string >,
          const std::string >(
          &tom::bodyOccultationViabilitySettings ),
          py::arg("link_end_id" ),
          py::arg("occulting_body" ),
          get_docstring("body_occultation_viability_list").c_str() );


    py::class_<tss::ObservationSimulationSettings<double>,
            std::shared_ptr<tss::ObservationSimulationSettings<double>>>(m, "ObservationSimulationSettings",
                                                                         get_docstring("ObservationSimulationSettings").c_str() )
            .def_property("noise_function",
                          &tss::ObservationSimulationSettings<double>::getObservationNoiseFunction,
                          py::overload_cast< const std::function< Eigen::VectorXd( const double ) >& >(
                                  &tss::ObservationSimulationSettings<double>::setObservationNoiseFunction ),
                          get_docstring("ObservationSimulationSettings.noise_function").c_str() );

    py::class_<tss::TabulatedObservationSimulationSettings<double>,
               std::shared_ptr<tss::TabulatedObservationSimulationSettings<double>>,
               tss::ObservationSimulationSettings<double> >(m, "TabulatedObservationSimulationSettings",
                                                            get_docstring("TabulatedObservationSimulationSettings").c_str() );


    m.def("tabulated_simulation_settings",
              &tss::tabulatedObservationSimulationSettings< double >,
          py::arg("observable_type"),
          py::arg("link_ends" ),
          py::arg("simulation_times" ),
          py::arg("reference_link_end_type" ) = tom::receiver,
          py::arg("viability_settings" ) = std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >( ),
          py::arg("noise_function" ) = nullptr,
          get_docstring("tabulated_simulation_settings").c_str() );


    m.def("tabulated_simulation_settings_list",
              &tss::createTabulatedObservationSimulationSettingsList< double >,
          py::arg("link_ends_per_observable"),
          py::arg("simulation_times" ),
          py::arg("reference_link_end_type" ) = tom::receiver,
          py::arg("viability_settings" ) = std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >( ),
          get_docstring("tabulated_simulation_settings_list").c_str() );

    m.def("add_gaussian_noise_to_settings",
          py::overload_cast<
          const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< double > > >&,
          const double,
          const tom::ObservableType >(
                &tss::addGaussianNoiseFunctionToObservationSimulationSettingsPy ),
            py::arg("observation_simulation_settings"),
            py::arg("noise_amplitude"),
            py::arg("observable_type"),
            get_docstring("add_gaussian_noise_to_settings").c_str() );


    m.def("add_viability_check_to_settings",
          py::overload_cast<
          const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< double > > >&,
          const std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >& >(
                &tss::addViabilityToObservationSimulationSettingsPy ),
            py::arg("observation_simulation_settings"),
            py::arg("viability_settings"),
            get_docstring("add_viability_check_to_settings", 0).c_str() );


    m.def("add_viability_check_to_settings",
          py::overload_cast<
          const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< double > > >&,
          const std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >&,
          const tom::ObservableType >(
                &tss::addViabilityToObservationSimulationSettingsPy ),
            py::arg("observation_simulation_settings"),
            py::arg("viability_settings"),
            py::arg("observable_type"),
            get_docstring("add_viability_check_to_settings", 1).c_str() );


    py::class_<tss::ObservationDependentVariableSettings,
            std::shared_ptr<tss::ObservationDependentVariableSettings>>(
                m, "ObservationDependentVariableSettings",
                get_docstring("ObservationDependentVariableSettings").c_str() );

    m.def("add_dependent_variables_to_settings",
          py::overload_cast<
          const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< double > > >&,
          const std::vector< std::shared_ptr< tss::ObservationDependentVariableSettings > >&,
          const tss::SystemOfBodies& >(
                &tss::addDependentVariablesToObservationSimulationSettings< double > ),
            py::arg("observation_simulation_settings"),
            py::arg("dependent_variable_settings" ),
            py::arg("bodies" ),
            get_docstring("add_dependent_variables_to_settings").c_str() );



//    m.def("gaussian_noise_function",
//              &ts::getGaussianDistributionNoiseFunction,
//          py::arg("standard_deviation"),
//          py::arg("mean") = 0.0,
//          py::arg("seed") = time(NULL),
//          py::arg("observable_size") = 1);
}

}// namespace observation
}// namespace estimation_setup
}// namespace numerical_simulation
}// namespace tudatpy

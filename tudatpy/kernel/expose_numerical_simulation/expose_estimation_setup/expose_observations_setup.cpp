/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_observations_setup.h"

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

void expose_observations_setup(py::module &m) {

    py::enum_< tom::LinkEndType >(m, "LinkEndType")
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
          py::arg("receivers") );

    m.def("one_way_uplink_link_ends",
          &tom::getOneWayUplinkLinkEndsList,
          py::arg("transmitters"),
          py::arg("receiver") );

    py::enum_< tom::ObservableType >(m, "ObservableType")
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
                m, "DopplerProperTimeRateSettings");

    py::class_<tom::ObservationModelSettings,
            std::shared_ptr<tom::ObservationModelSettings>>(
                m, "ObservationSettings");

    m.def("one_way_range",
          &tom::oneWayRangeSettings,
          py::arg("link_ends"),
          py::arg("light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
          py::arg("bias_settings") = nullptr );

    m.def("angular_position",
          &tom::angularPositionSettings,
          py::arg("link_ends"),
          py::arg("light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
          py::arg("bias_settings") = nullptr );

    m.def("cartesian_position",
          &tom::positionObservableSettings,
          py::arg("link_ends"),
          py::arg("bias_settings") = nullptr );

    m.def("cartesian_velocity",
          &tom::velocityObservableSettings,
          py::arg("link_ends"),
          py::arg("bias_settings") = nullptr );

    m.def("313_euler_angles",
          &tom::eulerAngle313ObservableSettings,
          py::arg("link_ends"),
          py::arg("bias_settings") = nullptr );

    m.def("one_way_open_loop_doppler",
          &tom::oneWayOpenLoopDoppler,
          py::arg("link_ends"),
          py::arg("light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
          py::arg("bias_settings") = nullptr,
          py::arg("transmitter_proper_time_rate_settings") = nullptr,
          py::arg("receiver_proper_time_rate_settings") = nullptr );

    m.def("two_way_open_loop_doppler",
          &tom::twoWayOpenLoopDoppler,
          py::arg("uplink_doppler_settings" ),
          py::arg("downlink_doppler_settings" ),
          py::arg("bias_settings") = nullptr );

    m.def("one_way_closed_loop_doppler",
          py::overload_cast<
          const tom::LinkEnds&,
          const double,
          const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >,
          const std::shared_ptr< tom::ObservationBiasSettings > >( &tom::oneWayClosedLoopDoppler ),
          py::arg("link_ends" ),
          py::arg("integration_time" ),
          py::arg("light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
          py::arg("bias_settings") = nullptr );

    m.def("one_way_closed_loop_doppler",
          py::overload_cast<
          const tom::LinkEnds&,
          const std::function< double( const double ) >,
          const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >,
          const std::shared_ptr< tom::ObservationBiasSettings > >( &tom::oneWayClosedLoopDoppler ),
          py::arg("link_ends" ),
          py::arg("integration_time_function" ),
          py::arg("light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
          py::arg("bias_settings") = nullptr );

    m.def("n_way_range",
          &tom::nWayRange,
          py::arg("one_way_range_settings" ),
          py::arg("bias_settings" ) = nullptr,
          py::arg("retransmission_times_function" ) = nullptr );

    py::class_<tom::LightTimeCorrectionSettings,
            std::shared_ptr<tom::LightTimeCorrectionSettings>>(
                m, "LightTimeCorrectionSettings");

    m.def("first_order_relativistic_light_time_correction",
          &tom::firstOrderRelativisticLightTimeCorrectionSettings,
          py::arg("perturbing_bodies") );

    py::class_<tom::ObservationBiasSettings,
            std::shared_ptr<tom::ObservationBiasSettings>>(
                m, "ObservationBiasSettings");

    m.def("bias",
          &tom::constantAbsoluteBias,
          py::arg("bias_value") );

    m.def("relative_bias",
          &tom::constantRelativeBias,
          py::arg("bias_value") );

    m.def("arcwise_bias",
          py::overload_cast<
          const std::vector< double >&,
          const std::vector< Eigen::VectorXd >&,
          const tom::LinkEndType >( &tom::arcWiseAbsoluteBias ),
          py::arg("arc_start_times" ),
          py::arg("bias_values"),
          py::arg("time_link_end" ) );

    m.def("arcwise_bias",
          py::overload_cast<
          const std::map< double, Eigen::VectorXd >&,
          const tom::LinkEndType >( &tom::arcWiseAbsoluteBias ),
          py::arg("bias_values_per_start_time"),
          py::arg("time_link_end" ) );

    m.def("arcwise_relative_bias",
          py::overload_cast<
          const std::vector< double >&,
          const std::vector< Eigen::VectorXd >&,
          const tom::LinkEndType >( &tom::arcWiseRelativeBias ),
          py::arg("arc_start_times" ),
          py::arg("bias_values"),
          py::arg("time_link_end" ) );

    m.def("arcwise_relative_bias",
          py::overload_cast<
          const std::map< double, Eigen::VectorXd >&,
          const tom::LinkEndType >( &tom::arcWiseRelativeBias ),
          py::arg("bias_values_per_start_time"),
          py::arg("time_link_end" ) );

    m.def("combined_bias",
          &tom::multipleObservationBiasSettings,
          py::arg("bias_list") );


    py::enum_< tom::ObservationViabilityType >(m, "ObservationViabilityType")
            .value("minimum_elevation_angle", tom::ObservationViabilityType::minimum_elevation_angle )
            .value("body_avoidance_angle", tom::ObservationViabilityType::body_avoidance_angle )
            .value("body_occultation", tom::ObservationViabilityType::body_occultation )
            .export_values();


    py::class_<tom::ObservationViabilitySettings,
            std::shared_ptr<tom::ObservationViabilitySettings>>(
                m, "ObservationViabilitySettings")
            .def(py::init< const tom::ObservationViabilityType,
                 const std::pair< std::string, std::string >,
                 const std::string,
                 const double >(),
                 py::arg("viability_type"),
                 py::arg("associated_link_end"),
                 py::arg("string_input"),
                 py::arg("double_input") );

    m.def("elevation_angle_viability_list",
          py::overload_cast<
          const std::vector< std::pair< std::string, std::string > >,
          const double >(
          &tom::elevationAngleViabilitySettings ),
          py::arg("link_ends_list" ),
          py::arg("elevation_angle" ) );

    m.def("elevation_angle_viability",
          py::overload_cast<
          const std::pair< std::string, std::string >,
          const double >(
          &tom::elevationAngleViabilitySettings ),
          py::arg("link_end" ),
          py::arg("elevation_angle" ) );


    m.def("body_avoidance_viability",
          py::overload_cast<
          const std::pair< std::string, std::string >,
          const std::string,
          const double >(
          &tom::bodyAvoidanceAngleViabilitySettings ),
          py::arg("link_end" ),
          py::arg("body_to_avoid" ),
          py::arg("avoidance_angle") );

    m.def("body_avoidance_viability_list",
          py::overload_cast<
          const std::vector< std::pair< std::string, std::string > >,
          const std::string,
          const double >(
          &tom::bodyAvoidanceAngleViabilitySettings ),
          py::arg("link_ends_list" ),
          py::arg("body_to_avoid" ),
          py::arg("avoidance_angle") );


    m.def("body_occultation_viability",
          py::overload_cast<
          const std::vector< std::pair< std::string, std::string > >,
          const std::string >(
          &tom::bodyOccultationViabilitySettings ),
          py::arg("link_ends_list" ),
          py::arg("occulting_body" ) );


    m.def("body_occultation_viability",
          py::overload_cast<
          const std::pair< std::string, std::string >,
          const std::string >(
          &tom::bodyOccultationViabilitySettings ),
          py::arg("link_end" ),
          py::arg("occulting_body" ) );

    m.def("create_observation_simulators",
          py::overload_cast< const std::vector< std::shared_ptr< tom::ObservationModelSettings > >&, const tss::SystemOfBodies& >(
              &tom::createObservationSimulators< double, double > ),
          py::arg( "observation_settings" ),
          py::arg( "bodies" ) );

    py::class_<tss::ObservationSimulationSettings<double>,
               std::shared_ptr<tss::ObservationSimulationSettings<double>>>(m, "ObservationSimulationSettings");

    py::class_<tss::TabulatedObservationSimulationSettings<double>,
               std::shared_ptr<tss::TabulatedObservationSimulationSettings<double>>,
               tss::ObservationSimulationSettings<double> >(m, "TabulatedObservationSimulationSettings")
            .def(py::init<
                 const tom::ObservableType, const tom::LinkEnds, const std::vector< double >, const tom::LinkEndType,
                 const std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >&,
                 const std::function< Eigen::VectorXd( const double ) > >(),
                 py::arg("observable_type"),
                 py::arg("link_ends"),
                 py::arg("observation_times"),
                 py::arg("reference_link_end") = tom::receiver,
                 py::arg("viability_settings") = std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >( ),
                 py::arg("noise_function") = nullptr );


    m.def("tabulated_simulation_settings",
              &tss::tabulatedObservationSimulationSettings< double >,
          py::arg("observable_type"),
          py::arg("link_ends" ),
          py::arg("simulation_times" ),
          py::arg("reference_link_end" ) = tom::receiver,
          py::arg("viability_settings" ) = std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >( ),
          py::arg("noise_function" ) = nullptr );


    m.def("create_tabulated_simulation_settings",
              &tss::createTabulatedObservationSimulationSettingsList< double >,
          py::arg("link_ends_per_observable"),
          py::arg("simulation_times" ) );

//    m.def("add_noise_to_settings",
//          py::overload_cast<
//          const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< double > > >&,
//          const std::function< Eigen::VectorXd( const double ) >,
//          const tom::ObservableType >(
//                &tss::addNoiseFunctionToObservationSimulationSettings< double, Eigen::VectorXd, const tom::ObservableType > ),
//            py::arg("observation_simulation_settings"),
//            py::arg("noise_functiton"),
//            py::arg("observable_type") );

    m.def("add_gaussian_noise_to_settings",
          py::overload_cast<
          const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< double > > >&,
          const double,
          const tom::ObservableType >(
                &tss::addGaussianNoiseFunctionToObservationSimulationSettingsPy ),
            py::arg("observation_simulation_settings"),
            py::arg("noise_amplitude"),
            py::arg("observable_type") );


    m.def("add_viability_check_to_settings",
          py::overload_cast<
          const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< double > > >&,
          const std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >& >(
                &tss::addViabilityToObservationSimulationSettingsPy ),
            py::arg("observation_simulation_settings"),
            py::arg("viability_settings") );

    m.def("add_viability_check_to_settings",
          py::overload_cast<
          const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< double > > >&,
          const std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >&,
          const tom::ObservableType >(
                &tss::addViabilityToObservationSimulationSettingsPy ),
            py::arg("observation_simulation_settings"),
            py::arg("viability_settings"),
            py::arg("observable_type") );

    py::class_<tss::ObservationDependentVariableSettings,
            std::shared_ptr<tss::ObservationDependentVariableSettings>>(
                m, "ObservationDependentVariableSettings");

//    m.def("add_dependent_variables_to_settings",
//          py::overload_cast<
//          const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< double > > >&,
//          const std::vector< std::shared_ptr< tss::ObservationDependentVariableSettings > >&,
//          const tss::SystemOfBodies& >(
//                &tss::addDependentVariablesToObservationSimulationSettings< double > ),
//            py::arg("observation_simulation_settings"),
//            py::arg("dependent_variable_settings" ),
//            py::arg("bodies" ) );


    m.def("simulate_observations",
              &tss::simulateObservations< >,
          py::arg("observation_to_simulate"),
          py::arg("observation_simulators" ),
          py::arg("bodies") );

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

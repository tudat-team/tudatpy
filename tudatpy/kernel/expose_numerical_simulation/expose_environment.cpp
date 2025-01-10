/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_environment.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/astro/aerodynamics.h>
#include <tudat/astro/ephemerides.h>
#include <tudat/astro/gravitation.h>
#include <tudat/basics/deprecationWarnings.h>

#include "docstrings.h"
#include "scalarTypes.h"
#include "tudat/astro/ground_stations/groundStation.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"

// namespace py = pybind11;
// namespace tba = tudat::basic_astrodynamics;
// namespace tss = tudat::simulation_setup;
// namespace tp = tudat::propagators;
// namespace tinterp = tudat::interpolators;
// namespace te = tudat::ephemerides;
// namespace tni = tudat::numerical_integrators;
// namespace trf = tudat::reference_frames;
// namespace tmrf = tudat::root_finders;

namespace py = pybind11;

namespace tba = tudat::basic_astrodynamics;
namespace ta = tudat::aerodynamics;
namespace tr = tudat::reference_frames;
namespace te = tudat::ephemerides;
namespace tgs = tudat::ground_stations;
namespace tr = tudat::reference_frames;
namespace tg = tudat::gravitation;
namespace trf = tudat::reference_frames;
namespace tss = tudat::simulation_setup;
namespace ti = tudat::interpolators;
namespace tsm = tudat::system_models;
namespace tom = tudat::observation_models;

namespace tudat
{

namespace aerodynamics
{

double getTotalSurfaceArea( const std::shared_ptr< HypersonicLocalInclinationAnalysis > coefficientGenerator )
{
    double totalSurfaceArea = 0.0;
    for( int i = 0; i < coefficientGenerator->getNumberOfVehicleParts( ); i++ )
    {
        totalSurfaceArea += std::fabs( coefficientGenerator->getVehiclePart( i )->getTotalArea( ) );
    }
    return totalSurfaceArea;
}

//! Function that saves the vehicle mesh data used for a HypersonicLocalInclinationAnalysis to a file
std::pair< std::vector< Eigen::Vector3d >, std::vector< Eigen::Vector3d > > getVehicleMesh(
        const std::shared_ptr< HypersonicLocalInclinationAnalysis > localInclinationAnalysis )
{
    std::vector< boost::multi_array< Eigen::Vector3d, 2 > > meshPoints = localInclinationAnalysis->getMeshPoints( );
    std::vector< boost::multi_array< Eigen::Vector3d, 2 > > meshSurfaceNormals = localInclinationAnalysis->getPanelSurfaceNormals( );

    //    boost::array< int, 3 > independentVariables;
    //    independentVariables[ 0 ] = 0;
    //    independentVariables[ 1 ] = 6;
    //    independentVariables[ 2 ] = 0;

    //    std::vector< std::vector< std::vector< double > > > pressureCoefficients =
    //            localInclinationAnalysis->getPressureCoefficientList( independentVariables );

    int counter = 0;
    std::vector< Eigen::Vector3d > meshPointsList;
    std::vector< Eigen::Vector3d > surfaceNormalsList;
    //    std::map< int, Eigen::Vector1d > pressureCoefficientsList;

    for( unsigned int i = 0; i < meshPoints.size( ); i++ )
    {
        for( unsigned int j = 0; j < meshPoints.at( i ).shape( )[ 0 ] - 1; j++ )
        {
            for( unsigned int k = 0; k < meshPoints.at( i ).shape( )[ 1 ] - 1; k++ )
            {
                meshPointsList.push_back( meshPoints[ i ][ j ][ k ] );
                surfaceNormalsList.push_back( meshSurfaceNormals[ i ][ j ][ k ] );
                //                pressureCoefficientsList[ counter ] = ( Eigen::Vector1d( ) << pressureCoefficients[ i ][ j ][ k ]
                //                ).finished( );
                counter++;
            }
        }
    }

    return std::make_pair( meshPointsList, surfaceNormalsList );
}

}  // namespace aerodynamics

}  // namespace tudat

namespace tudatpy
{
namespace numerical_simulation
{
namespace environment
{

void expose_environment( py::module &m )
{
    py::enum_< ta::AerodynamicCoefficientsIndependentVariables >(
            m, "AerodynamicCoefficientsIndependentVariables", get_docstring( "AerodynamicCoefficientsIndependentVariables" ).c_str( ) )
            .value( "mach_number_dependent",
                    ta::AerodynamicCoefficientsIndependentVariables::mach_number_dependent,
                    get_docstring( "AerodynamicCoefficientsIndependentVariables.mach_number_dependent" ).c_str( ) )
            .value( "angle_of_attack_dependent",
                    ta::AerodynamicCoefficientsIndependentVariables::angle_of_attack_dependent,
                    get_docstring( "AerodynamicCoefficientsIndependentVariables.angle_of_attack_dependent" ).c_str( ) )
            .value( "sideslip_angle_dependent",
                    ta::AerodynamicCoefficientsIndependentVariables::angle_of_sideslip_dependent,
                    get_docstring( "AerodynamicCoefficientsIndependentVariables.sideslip_angle_dependent" ).c_str( ) )
            .value( "altitude_dependent",
                    ta::AerodynamicCoefficientsIndependentVariables::altitude_dependent,
                    get_docstring( "AerodynamicCoefficientsIndependentVariables.altitude_dependent" ).c_str( ) )
            .value( "time_dependent",
                    ta::AerodynamicCoefficientsIndependentVariables::time_dependent,
                    get_docstring( "AerodynamicCoefficientsIndependentVariables.time_dependent" ).c_str( ) )
            .value( "temperature_dependent",
                    ta::AerodynamicCoefficientsIndependentVariables::temperature_dependent,
                    get_docstring( "AerodynamicCoefficientsIndependentVariables.temperature_dependent" ).c_str( ) )
            .value( "velocity_dependent",
                    ta::AerodynamicCoefficientsIndependentVariables::velocity_dependent,
                    get_docstring( "AerodynamicCoefficientsIndependentVariables.velocity_dependent" ).c_str( ) )
            .value( "he_number_density_dependent",
                    ta::AerodynamicCoefficientsIndependentVariables::he_number_density_dependent,
                    get_docstring( "AerodynamicCoefficientsIndependentVariables.he_number_density_dependent" ).c_str( ) )
            .value( "o_number_density_dependent",
                    ta::AerodynamicCoefficientsIndependentVariables::o_number_density_dependent,
                    get_docstring( "AerodynamicCoefficientsIndependentVariables.o_number_density_dependent" ).c_str( ) )
            .value( "n2_number_density_dependent",
                    ta::AerodynamicCoefficientsIndependentVariables::n2_number_density_dependent,
                    get_docstring( "AerodynamicCoefficientsIndependentVariables.n2_number_density_dependent" ).c_str( ) )
            .value( "o2_number_density_dependent",
                    ta::AerodynamicCoefficientsIndependentVariables::o2_number_density_dependent,
                    get_docstring( "AerodynamicCoefficientsIndependentVariables.o2_number_density_dependent" ).c_str( ) )
            .value( "ar_number_density_dependent",
                    ta::AerodynamicCoefficientsIndependentVariables::ar_number_density_dependent,
                    get_docstring( "AerodynamicCoefficientsIndependentVariables.ar_number_density_dependent" ).c_str( ) )
            .value( "h_number_density_dependent",
                    ta::AerodynamicCoefficientsIndependentVariables::h_number_density_dependent,
                    get_docstring( "AerodynamicCoefficientsIndependentVariables.h_number_density_dependent" ).c_str( ) )
            .value( "n_number_density_dependent",
                    ta::AerodynamicCoefficientsIndependentVariables::n_number_density_dependent,
                    get_docstring( "AerodynamicCoefficientsIndependentVariables.n_number_density_dependent" ).c_str( ) )
            .value( "anomalous_o_number_density_dependent",
                    ta::AerodynamicCoefficientsIndependentVariables::anomalous_o_number_density_dependent,
                    get_docstring( "AerodynamicCoefficientsIndependentVariables.anomalous_o_number_density_dependent" ).c_str( ) )
            .value( "control_surface_deflection_dependent",
                    ta::AerodynamicCoefficientsIndependentVariables::control_surface_deflection_dependent,
                    get_docstring( "AerodynamicCoefficientsIndependentVariables.control_surface_deflection_dependent" ).c_str( ) )
            .value( "undefined_independent_variable",
                    ta::AerodynamicCoefficientsIndependentVariables::undefined_independent_variable,
                    get_docstring( "AerodynamicCoefficientsIndependentVariables.undefined_independent_variable" ).c_str( ) )
            .export_values( );

    py::enum_< ta::AerodynamicCoefficientFrames >(
            m, "AerodynamicCoefficientFrames", get_docstring( "AerodynamicCoefficientFrames" ).c_str( ) )
            .value( "positive_body_fixed_frame_coefficients",
                    ta::AerodynamicCoefficientFrames::body_fixed_frame_coefficients,
                    get_docstring( "AerodynamicCoefficientFrames.positive_body_fixed_frame_coefficients" ).c_str( ) )
            .value( "negative_body_fixed_frame_coefficients",
                    ta::AerodynamicCoefficientFrames::negative_body_fixed_frame_coefficients,
                    get_docstring( "AerodynamicCoefficientFrames.negative_body_fixed_frame_coefficients" ).c_str( ) )
            .value( "positive_aerodynamic_frame_coefficients",
                    ta::AerodynamicCoefficientFrames::positive_aerodynamic_frame_coefficients,
                    get_docstring( "AerodynamicCoefficientFrames.positive_aerodynamic_frame_coefficients" ).c_str( ) )
            .value( "negative_aerodynamic_frame_coefficients",
                    ta::AerodynamicCoefficientFrames::negative_aerodynamic_frame_coefficients,
                    get_docstring( "AerodynamicCoefficientFrames.negative_aerodynamic_frame_coefficients" ).c_str( ) )
            .export_values( );

    py::class_< ta::AtmosphereModel, std::shared_ptr< ta::AtmosphereModel > >(
            m, "AtmosphereModel", get_docstring( "AtmosphereModel" ).c_str( ) )
            .def( "get_density",
                  &ta::AtmosphereModel::getDensity,
                  py::arg( "altitude" ),
                  py::arg( "longitude" ),
                  py::arg( "latitude" ),
                  py::arg( "time" ),
                  get_docstring( "AtmosphereModel.get_density" ).c_str( ) )
            .def( "get_pressure",
                  &ta::AtmosphereModel::getPressure,
                  py::arg( "altitude" ),
                  py::arg( "longitude" ),
                  py::arg( "latitude" ),
                  py::arg( "time" ),
                  get_docstring( "AtmosphereModel.get_pressure" ).c_str( ) )
            .def( "get_temperature",
                  &ta::AtmosphereModel::getTemperature,
                  py::arg( "altitude" ),
                  py::arg( "longitude" ),
                  py::arg( "latitude" ),
                  py::arg( "time" ),
                  get_docstring( "AtmosphereModel.get_temperature" ).c_str( ) )
            .def( "get_speed_of_sound",
                  &ta::AtmosphereModel::getSpeedOfSound,
                  py::arg( "altitude" ),
                  py::arg( "longitude" ),
                  py::arg( "latitude" ),
                  py::arg( "time" ),
                  get_docstring( "AtmosphereModel.get_speed_of_sound" ).c_str( ) )
            .def( "get_number_density",
                  &ta::AtmosphereModel::getNumberDensity,
                  py::arg( "species" ),
                  py::arg( "altitude" ),
                  py::arg( "longitude" ),
                  py::arg( "latitude" ),
                  py::arg( "time" ),
                  get_docstring( "AtmosphereModel.get_number_density" ).c_str( ) );

    py::class_< ta::AerodynamicCoefficientInterface, std::shared_ptr< ta::AerodynamicCoefficientInterface > >(
            m, "AerodynamicCoefficientInterface", get_docstring( "AerodynamicCoefficientInterface" ).c_str( ) )
            .def_property_readonly( "reference_area",
                                    &ta::AerodynamicCoefficientInterface::getReferenceArea,
                                    get_docstring( "AerodynamicCoefficientInterface.reference_area" ).c_str( ) )
            .def_property_readonly( "current_force_coefficients",
                                    &ta::AerodynamicCoefficientInterface::getCurrentForceCoefficients,
                                    get_docstring( "AerodynamicCoefficientInterface.current_force_coefficients" ).c_str( ) )
            .def_property_readonly( "current_moment_coefficients",
                                    &ta::AerodynamicCoefficientInterface::getCurrentMomentCoefficients,
                                    get_docstring( "AerodynamicCoefficientInterface.current_moment_coefficients" ).c_str( ) )
            .def_property_readonly( "current_coefficients",
                                    &ta::AerodynamicCoefficientInterface::getCurrentAerodynamicCoefficients,
                                    get_docstring( "AerodynamicCoefficientInterface.current_coefficients" ).c_str( ) )
            .def_property_readonly( "force_coefficient_frame",
                                    &ta::AerodynamicCoefficientInterface::getForceCoefficientsFrame,
                                    get_docstring( "AerodynamicCoefficientInterface.force_coefficient_frame" ).c_str( ) )
            .def_property_readonly( "moment_coefficient_frame",
                                    &ta::AerodynamicCoefficientInterface::getMomentCoefficientsFrame,
                                    get_docstring( "AerodynamicCoefficientInterface.moment_coefficient_frame" ).c_str( ) )
            .def_property_readonly( "independent_variable_names",
                                    &ta::AerodynamicCoefficientInterface::getIndependentVariableNames,
                                    get_docstring( "AerodynamicCoefficientInterface.independent_variable_names" ).c_str( ) )
            .def_property_readonly(
                    "current_control_surface_free_force_coefficients",
                    &ta::AerodynamicCoefficientInterface::getCurrentControlSurfaceFreeForceCoefficients,
                    get_docstring( "AerodynamicCoefficientInterface.current_control_surface_free_force_coefficients" ).c_str( ) )
            .def_property_readonly(
                    "current_control_surface_free_moment_coefficients",
                    &ta::AerodynamicCoefficientInterface::getCurrentControlSurfaceFreeMomentCoefficients,
                    get_docstring( "AerodynamicCoefficientInterface.current_control_surface_free_moment_coefficients" ).c_str( ) )
            .def_property_readonly( "control_surface_independent_variable_names",
                                    &ta::AerodynamicCoefficientInterface::getControlSurfaceIndependentVariables,
                                    get_docstring( "AerodynamicCoefficientInterface.control_surface_independent_variable_names" ).c_str( ) )
            .def( "current_control_surface_force_coefficient_increment",
                  &ta::AerodynamicCoefficientInterface::getCurrentForceCoefficientIncrement,
                  py::arg( "control_surface_name" ),
                  get_docstring( "AerodynamicCoefficientInterface.current_control_surface_force_coefficient_increment" ).c_str( ) )
            .def( "current_control_surface_moment_coefficient_increment",
                  &ta::AerodynamicCoefficientInterface::getCurrentMomentCoefficientIncrement,
                  py::arg( "control_surface_name" ),
                  get_docstring( "AerodynamicCoefficientInterface.current_control_surface_moment_coefficient_increment" ).c_str( ) )
            .def( "set_control_surface_increments",
                  &ta::AerodynamicCoefficientInterface::setControlSurfaceIncrements,
                  py::arg( "control_surface_list" ),
                  get_docstring( "AerodynamicCoefficientInterface.set_control_surface_increments" ).c_str( ) )
            .def( "update_coefficients",
                  &ta::AerodynamicCoefficientInterface::updateCurrentCoefficients,
                  py::arg( "independent_variables" ),
                  py::arg( "time" ),
                  get_docstring( "AerodynamicCoefficientInterface.update_coefficients" ).c_str( ) )
            .def( "update_full_coefficients",
                  &ta::AerodynamicCoefficientInterface::updateFullCurrentCoefficients,
                  py::arg( "independent_variables" ),
                  py::arg( "control_surface_independent_variables" ),
                  py::arg( "time" ),
                  py::arg( "check_force_contribution" ) = true,
                  get_docstring( "AerodynamicCoefficientInterface.update_full_coefficients" ).c_str( ) );

    py::class_< ta::AerodynamicCoefficientGenerator< 3, 6 >,
                std::shared_ptr< ta::AerodynamicCoefficientGenerator< 3, 6 > >,
                ta::AerodynamicCoefficientInterface >( m, "AerodynamicCoefficientGenerator36", "<no_doc, only_dec>" );

    py::class_< ta::HypersonicLocalInclinationAnalysis,
                std::shared_ptr< ta::HypersonicLocalInclinationAnalysis >,
                ta::AerodynamicCoefficientGenerator< 3, 6 > >( m, "HypersonicLocalInclinationAnalysis" )
            .def( py::init< const std::vector< std::vector< double > > &,
                            const std::shared_ptr< tudat::SurfaceGeometry >,
                            const std::vector< int > &,
                            const std::vector< int > &,
                            const std::vector< bool > &,
                            const std::vector< std::vector< int > > &,
                            const double,
                            const double,
                            const Eigen::Vector3d &,
                            const bool >( ),
                  py::arg( "independent_variable_points" ),
                  py::arg( "body_shape" ),
                  py::arg( "number_of_lines" ),
                  py::arg( "number_of_points" ),
                  py::arg( "invert_orders" ),
                  py::arg( "selected_methods" ),
                  py::arg( "reference_area" ),
                  py::arg( "reference_length" ),
                  py::arg( "moment_reference_point" ),
                  py::arg( "save_pressure_coefficients" ) = false,
                  get_docstring( "HypersonicLocalInclinationAnalysis.ctor" ).c_str( ) )
            .def( "clear_data", &ta::HypersonicLocalInclinationAnalysis::clearData );

    py::class_< ta::ControlSurfaceIncrementAerodynamicInterface, std::shared_ptr< ta::ControlSurfaceIncrementAerodynamicInterface > >(
            m, "ControlSurfaceIncrementAerodynamicInterface", "<no_doc, only_dec>" );

    py::class_< ta::CustomControlSurfaceIncrementAerodynamicInterface,
                std::shared_ptr< ta::CustomControlSurfaceIncrementAerodynamicInterface >,
                ta::ControlSurfaceIncrementAerodynamicInterface >(
            m, "CustomControlSurfaceIncrementAerodynamicInterface", "<no_doc, only_dec>" )
            .def( py::init< const std::function< Eigen::Vector6d( const std::vector< double > & ) >,
                            const std::vector< ta::AerodynamicCoefficientsIndependentVariables > >( ),
                  py::arg( "coefficient_function" ),
                  py::arg( "independent_variable_names" ) );

    m.def( "get_default_local_inclination_mach_points",
           &ta::getDefaultHypersonicLocalInclinationMachPoints,
           py::arg( "mach_regime" ) = "Full" );

    m.def( "get_default_local_inclination_angle_of_attack_points", &ta::getDefaultHypersonicLocalInclinationAngleOfAttackPoints );

    m.def( "get_default_local_inclination_sideslip_angle_points", &ta::getDefaultHypersonicLocalInclinationAngleOfSideslipPoints );

    m.def( "save_vehicle_mesh_to_file",
           &ta::saveVehicleMeshToFile,
           py::arg( "local_inclination_analysis_object" ),
           py::arg( "output_directory" ),
           py::arg( "output_file_prefix" ) = "",
           get_docstring( "save_vehicle_mesh_to_file" ).c_str( ) );

    m.def( "get_local_inclination_total_vehicle_area", &ta::getTotalSurfaceArea, py::arg( "local_inclination_analysis_object" ) );

    m.def( "get_local_inclination_mesh", &ta::getVehicleMesh, py::arg( "local_inclination_analysis_object" ) );

    py::class_< tsm::VehicleSystems, std::shared_ptr< tsm::VehicleSystems > >(
            m, "VehicleSystems", get_docstring( "VehicleSystems" ).c_str( ) )
            .def( py::init<>( ) )
            .def( "set_control_surface_deflection",
                  &tsm::VehicleSystems::setCurrentControlSurfaceDeflection,
                  py::arg( "control_surface_id" ),
                  py::arg( "deflection_angle" ),
                  get_docstring( "VehicleSystems.set_control_surface_deflection" ).c_str( ) )
            .def( "set_transponder_turnaround_ratio",
                  py::overload_cast< std::map< std::pair< tom::FrequencyBands, tom::FrequencyBands >, double > & >(
                          &tsm::VehicleSystems::setTransponderTurnaroundRatio ),
                  py::arg( "transponder_ratio_per_uplink_and_downlink_frequency_band" ),
                  get_docstring( "VehicleSystems.set_transponder_turnaround_ratio" ).c_str( ) )
            .def( "set_default_transponder_turnaround_ratio_function",
                  &tsm::VehicleSystems::setDefaultTransponderTurnaroundRatio,
                  get_docstring( "VehicleSystems.set_default_transponder_turnaround_ratio_function" ).c_str( ) )
            .def( "get_control_surface_deflection",
                  &tsm::VehicleSystems::getCurrentControlSurfaceDeflection,
                  py::arg( "control_surface_id" ),
                  get_docstring( "VehicleSystems.get_control_surface_deflection" ).c_str( ) )
            .def( "set_reference_point",
                  py::overload_cast< const std::string, const Eigen::Vector3d &, const std::string, const std::string >(
                          &tsm::VehicleSystems::setReferencePointPosition ),
                  py::arg( "reference_point" ),
                  py::arg( "location" ),
                  py::arg( "frame_origin" ) = "",
                  py::arg( "frame_orientation" ) = "",
                  get_docstring( "VehicleSystems.set_reference_point" ).c_str( ) )
            .def( "set_reference_point",
                  py::overload_cast< const std::string, std::shared_ptr< te::Ephemeris > >(
                          &tsm::VehicleSystems::setReferencePointPosition ),
                  py::arg( "reference_point" ),
                  py::arg( "ephemeris" ),
                  get_docstring( "VehicleSystems.set_reference_point" ).c_str( ) )

            .def( "get_engine_model",
                  &tsm::VehicleSystems::getEngineModel,
                  py::arg( "engine_name" ),
                  get_docstring( "VehicleSystems.get_engine_model" ).c_str( ) )
            .def( "set_timing_system", &tsm::VehicleSystems::setTimingSystem, py::arg( "timing_system" ) );

    py::class_< tsm::TimingSystem, std::shared_ptr< tsm::TimingSystem > >( m,
                                                                           "TimingSystem",
                                                                           get_docstring( "TimingSystem" ).c_str( ) )

            .def(  // ctor 1
                    py::init< const std::vector< tudat::Time >,
                              const std::vector< double >,
                              const std::function< std::function< double( const double ) >( const double, const double, const double ) >,
                              const double >( ),
                    py::arg( "arc_times" ),
                    py::arg( "all_arcs_polynomial_drift_coefficients" ) = std::vector< double >( ),
                    py::arg( "clock_noise_generation_function" ) = nullptr,
                    py::arg( "clock_noise_time_step" ) = 1.0E-3 )
            .def(  // ctor 2
                    py::init< const std::vector< tudat::Time >,
                              const std::vector< std::vector< double > >,
                              const std::function< std::function< double( const double ) >( const double, const double, const double ) >,
                              const double >( ),
                    py::arg( "arc_times" ),
                    py::arg( "polynomial_drift_coefficients" ),
                    py::arg( "clock_noise_generation_function" ) = nullptr,
                    py::arg( "clock_noise_time_step" ) = 1.0E-3 )
            .def(  // ctor 3
                    py::init< const std::vector< std::vector< double > >,
                              const std::vector< std::function< double( const double ) > >,
                              const std::vector< tudat::Time > >( ),
                    py::arg( "polynomial_drift_coefficients" ),
                    py::arg( "stochastic_clock_noise_functions" ),
                    py::arg( "arc_times" ) );

    py::class_< tsm::EngineModel, std::shared_ptr< tsm::EngineModel > >( m, "EngineModel" )
            .def_property_readonly( "thrust_magnitude_calculator", &tsm::EngineModel::getThrustMagnitudeWrapper );

    /*!
     **************   FLIGHT CONDITIONS AND ASSOCIATED FUNCTIONALITY  ******************
     */

    py::enum_< trf::AerodynamicsReferenceFrameAngles >( m, "AerodynamicsReferenceFrameAngles" )
            .value( "latitude_angle", trf::AerodynamicsReferenceFrameAngles::latitude_angle )
            .value( "longitude_angle", trf::AerodynamicsReferenceFrameAngles::longitude_angle )
            .value( "heading_angle", trf::AerodynamicsReferenceFrameAngles::heading_angle )
            .value( "flight_path_angle", trf::AerodynamicsReferenceFrameAngles::flight_path_angle )
            .value( "angle_of_attack", trf::AerodynamicsReferenceFrameAngles::angle_of_attack )
            .value( "angle_of_sideslip", trf::AerodynamicsReferenceFrameAngles::angle_of_sideslip )
            .value( "bank_angle", trf::AerodynamicsReferenceFrameAngles::bank_angle )
            .export_values( );

    py::enum_< trf::AerodynamicsReferenceFrames >(
            m, "AerodynamicsReferenceFrames", get_docstring( "AerodynamicsReferenceFrames" ).c_str( ) )
            .value( "inertial_frame",
                    trf::AerodynamicsReferenceFrames::inertial_frame,
                    get_docstring( "AerodynamicsReferenceFrames.inertial_frame" ).c_str( ) )
            .value( "corotating_frame",
                    trf::AerodynamicsReferenceFrames::corotating_frame,
                    get_docstring( "AerodynamicsReferenceFrames.corotating_frame" ).c_str( ) )
            .value( "vertical_frame",
                    trf::AerodynamicsReferenceFrames::vertical_frame,
                    get_docstring( "AerodynamicsReferenceFrames.vertical_frame" ).c_str( ) )
            .value( "trajectory_frame",
                    trf::AerodynamicsReferenceFrames::trajectory_frame,
                    get_docstring( "AerodynamicsReferenceFrames.trajectory_frame" ).c_str( ) )
            .value( "aerodynamic_frame",
                    trf::AerodynamicsReferenceFrames::aerodynamic_frame,
                    get_docstring( "AerodynamicsReferenceFrames.aerodynamic_frame" ).c_str( ) )
            .value( "body_frame",
                    trf::AerodynamicsReferenceFrames::body_frame,
                    get_docstring( "AerodynamicsReferenceFrames.body_frame" ).c_str( ) )
            .export_values( );

    py::class_< trf::AerodynamicAngleCalculator, std::shared_ptr< trf::AerodynamicAngleCalculator > >(
            m, "AerodynamicAngleCalculator", get_docstring( "AerodynamicAngleCalculator" ).c_str( ) )
            .def( "get_rotation_matrix_between_frames",
                  &trf::AerodynamicAngleCalculator::getRotationMatrixBetweenFrames,
                  py::arg( "original_frame" ),
                  py::arg( "target_frame" ),
                  get_docstring( "AerodynamicAngleCalculator.get_rotation_matrix_between_frames" ).c_str( ) )
            .def( "get_angle",
                  &trf::AerodynamicAngleCalculator::getAerodynamicAngle,
                  py::arg( "angle_type" ),
                  get_docstring( "AerodynamicAngleCalculator.get_angle" ).c_str( ) )
            // Function removed; error is shown
            .def( "set_body_orientation_angles",
                  &trf::AerodynamicAngleCalculator::setOrientationAngleFunctionsRemoved2,
                  py::arg( "angle_of_attack" ) = TUDAT_NAN,
                  py::arg( "angle_of_sideslip" ) = TUDAT_NAN,
                  py::arg( "bank_angle" ) = TUDAT_NAN,
                  py::arg( "silence_warnings" ) = false )
            // Function removed; error is shown
            .def( "set_body_orientation_angle_functions",
                  &trf::AerodynamicAngleCalculator::setOrientationAngleFunctionsRemoved1,
                  py::arg( "angle_of_attack_function" ) = std::function< double( ) >( ),    // <pybind11/functional.h>
                  py::arg( "angle_of_sideslip_function" ) = std::function< double( ) >( ),  // <pybind11/functional.h>
                  py::arg( "bank_angle_function" ) = std::function< double( ) >( ),         // <pybind11/functional.h>
                  py::arg( "angle_update_function" ) = std::function< void( const double ) >( ),
                  py::arg( "silence_warnings" ) = false );

    py::class_< ta::FlightConditions, std::shared_ptr< ta::FlightConditions > >(
            m, "FlightConditions", get_docstring( "FlightConditions" ).c_str( ) )
            //            .def(py::init<
            //                 const std::shared_ptr<tudat::basic_astrodynamics::BodyShapeModel>,
            //                 const std::shared_ptr<tudat::reference_frames::AerodynamicAngleCalculator>>(),
            //                 py::arg("shape_model"),
            //                 py::arg("aerodynamic_angle_calculator") = std::shared_ptr< tr::AerodynamicAngleCalculator>())
            .def( "update_conditions", &ta::FlightConditions::updateConditions, py::arg( "current_time" ) )
            .def_property_readonly( "aerodynamic_angle_calculator",
                                    &ta::FlightConditions::getAerodynamicAngleCalculator,
                                    get_docstring( "FlightConditions.aerodynamic_angle_calculator" ).c_str( ) )
            .def_property_readonly(
                    "longitude", &ta::FlightConditions::getCurrentLongitude, get_docstring( "FlightConditions.longitude" ).c_str( ) )
            .def_property_readonly(
                    "latitude", &ta::FlightConditions::getCurrentLatitude, get_docstring( "FlightConditions.latitude" ).c_str( ) )
            .def_property_readonly( "geodetic_latitude",
                                    &ta::FlightConditions::getCurrentGeodeticLatitude,
                                    get_docstring( "FlightConditions.latitude" ).c_str( ) )
            .def_property_readonly( "time", &ta::FlightConditions::getCurrentTime, get_docstring( "FlightConditions.time" ).c_str( ) )
            .def_property_readonly( "body_centered_body_fixed_state",
                                    &ta::FlightConditions::getCurrentBodyCenteredBodyFixedState,
                                    get_docstring( "FlightConditions.body_centered_body_fixed_state" ).c_str( ) )
            .def_property_readonly(
                    "altitude", &ta::FlightConditions::getCurrentAltitude, get_docstring( "FlightConditions.time" ).c_str( ) );

    py::class_< ta::AtmosphericFlightConditions, std::shared_ptr< ta::AtmosphericFlightConditions >, ta::FlightConditions >(
            m, "AtmosphericFlightConditions", get_docstring( "AtmosphericFlightConditions" ).c_str( ) )
            .def_property_readonly( "density",
                                    &ta::AtmosphericFlightConditions::getCurrentDensity,
                                    get_docstring( "AtmosphericFlightConditions.density" ).c_str( ) )
            .def_property_readonly( "temperature",
                                    &ta::AtmosphericFlightConditions::getCurrentFreestreamTemperature,
                                    get_docstring( "AtmosphericFlightConditions.temperature" ).c_str( ) )
            .def_property_readonly( "dynamic_pressure",
                                    &ta::AtmosphericFlightConditions::getCurrentDynamicPressure,
                                    get_docstring( "AtmosphericFlightConditions.dynamic_pressure" ).c_str( ) )
            .def_property_readonly( "pressure",
                                    &ta::AtmosphericFlightConditions::getCurrentPressure,
                                    get_docstring( "AtmosphericFlightConditions.pressure" ).c_str( ) )
            .def_property_readonly( "airspeed",
                                    &ta::AtmosphericFlightConditions::getCurrentAirspeed,
                                    get_docstring( "AtmosphericFlightConditions.airspeed" ).c_str( ) )
            .def_property_readonly( "mach_number",
                                    &ta::AtmosphericFlightConditions::getCurrentMachNumber,
                                    get_docstring( "AtmosphericFlightConditions.mach_number" ).c_str( ) )
            .def_property_readonly( "airspeed_velocity",
                                    &ta::AtmosphericFlightConditions::getCurrentAirspeedBasedVelocity,
                                    get_docstring( "AtmosphericFlightConditions.airspeed_velocity" ).c_str( ) )
            .def_property_readonly( "speed_of_sound",
                                    &ta::AtmosphericFlightConditions::getCurrentSpeedOfSound,
                                    get_docstring( "AtmosphericFlightConditions.speed_of_sound" ).c_str( ) )
            .def_property_readonly( "aero_coefficient_independent_variables",
                                    &ta::AtmosphericFlightConditions::getAerodynamicCoefficientIndependentVariables,
                                    get_docstring( "AtmosphericFlightConditions.aero_coefficient_independent_variables" ).c_str( ) )
            .def_property_readonly(
                    "control_surface_aero_coefficient_independent_variables",
                    &ta::AtmosphericFlightConditions::getControlSurfaceAerodynamicCoefficientIndependentVariables,
                    get_docstring( "AtmosphericFlightConditions.control_surface_aero_coefficient_independent_variables" ).c_str( ) )
            .def_property_readonly( "aerodynamic_coefficient_interface",
                                    &ta::AtmosphericFlightConditions::getAerodynamicCoefficientInterface,
                                    get_docstring( "AtmosphericFlightConditions.aerodynamic_coefficient_interface" ).c_str( ) );

    /*!
     **************   EPHEMERIDES  ******************
     */

    py::class_< te::Ephemeris, std::shared_ptr< te::Ephemeris > >( m, "Ephemeris", get_docstring( "Ephemeris" ).c_str( ) )
            .def( "cartesian_state",
                  &te::Ephemeris::getCartesianState,
                  py::arg( "current_time" ),
                  get_docstring( "Ephemeris.cartesian_state" ).c_str( ) )
            .def( "cartesian_position",
                  &te::Ephemeris::getCartesianPosition,
                  py::arg( "current_time" ),
                  get_docstring( "Ephemeris.cartesian_position" ).c_str( ) )
            .def( "cartesian_velocity",
                  &te::Ephemeris::getCartesianVelocity,
                  py::arg( "current_time" ),
                  get_docstring( "Ephemeris.cartesian_velocity" ).c_str( ) )
            .def_property_readonly(
                    "frame_origin", &te::Ephemeris::getReferenceFrameOrigin, get_docstring( "Ephemeris.frame_origin" ).c_str( ) )
            .def_property_readonly( "frame_orientation",
                                    &te::Ephemeris::getReferenceFrameOrientation,
                                    get_docstring( "Ephemeris.frame_orientation" ).c_str( ) );

    py::class_< te::ConstantEphemeris, std::shared_ptr< te::ConstantEphemeris >, te::Ephemeris >(
            m, "ConstantEphemeris", get_docstring( "ConstantEphemeris" ).c_str( ) )
            .def( py::init< const std::function< Eigen::Vector6d( ) >,  //<pybind11/functional.h>,<pybind11/eigen.h>
                            const std::string &,
                            const std::string & >( ),
                  py::arg( "constant_state_function" ),
                  py::arg( "reference_frame_origin" ) = "SSB",
                  py::arg( "reference_frame_orientation" ) = "ECLIPJ2000" )
            .def( py::init< const Eigen::Vector6d,  //<pybind11/eigen.h>
                            const std::string &,
                            const std::string & >( ),
                  py::arg( "constant_state" ),
                  py::arg( "reference_frame_origin" ) = "SSB",
                  py::arg( "reference_frame_orientation" ) = "ECLIPJ2000" )
            .def( "update_constant_state",
                  &te::ConstantEphemeris::updateConstantState,
                  py::arg( "new_state" ),
                  get_docstring( "ConstantEphemeris.update_constant_state" ).c_str( ) );

    py::class_< te::KeplerEphemeris, std::shared_ptr< te::KeplerEphemeris >, te::Ephemeris >( m, "KeplerEphemeris" );

    py::class_< te::MultiArcEphemeris, std::shared_ptr< te::MultiArcEphemeris >, te::Ephemeris >( m, "MultiArcEphemeris" )
            .def( py::init< const std::map< double, std::shared_ptr< te::Ephemeris > > &, const std::string &, const std::string & >( ),
                  py::arg( "single_arc_ephemerides" ),
                  py::arg( "reference_frame_origin" ) = "SSB",
                  py::arg( "reference_frame_orientation" ) = "ECLIPJ2000" );

    py::class_< te::TabulatedCartesianEphemeris< double, double >,
                std::shared_ptr< te::TabulatedCartesianEphemeris< double, double > >,
                te::Ephemeris >( m, "TabulatedEphemeris" )
            .def_property( "interpolator",
                           &te::TabulatedCartesianEphemeris< double, double >::getDynamicVectorSizeInterpolator,
                           py::overload_cast< const std::shared_ptr< ti::OneDimensionalInterpolator< double, Eigen::VectorXd > > >(
                                   &te::TabulatedCartesianEphemeris< double, double >::resetInterpolator ) );

    py::class_< te::Tle, std::shared_ptr< te::Tle > >( m, "Tle" )
            .def( py::init<  // ctor 1
                          const std::string & >( ),
                  py::arg( "lines" ) )
            .def( py::init<  // ctor 2
                          const std::string &,
                          const std::string & >( ),
                  py::arg( "line_1" ),
                  py::arg( "line_2" ) )
            .def( "get_epoch", &te::Tle::getEpoch )
            .def( "get_b_star", &te::Tle::getBStar )
            .def( "get_epoch", &te::Tle::getEpoch )
            .def( "get_inclination", &te::Tle::getInclination )
            .def( "get_right_ascension", &te::Tle::getRightAscension )
            .def( "get_eccentricity", &te::Tle::getEccentricity )
            .def( "get_arg_of_perigee", &te::Tle::getArgOfPerigee )
            .def( "get_mean_anomaly", &te::Tle::getMeanAnomaly )
            .def( "get_mean_motion", &te::Tle::getMeanMotion );

    py::class_< te::TleEphemeris, std::shared_ptr< te::TleEphemeris >, te::Ephemeris >( m, "TleEphemeris" )
            .def( py::init< const std::string &, const std::string &, const std::shared_ptr< te::Tle >, const bool >( ),
                  py::arg( "frame_origin" ) = "Earth",
                  py::arg( "frame_orientation" ) = "J2000",
                  py::arg( "tle" ) = nullptr,
                  py::arg( "use_sdp" ) = false );

    /*!
     **************   ROTATION MODELS  ******************
     */

    py::class_< te::RotationalEphemeris, std::shared_ptr< te::RotationalEphemeris > >(
            m, "RotationalEphemeris", get_docstring( "RotationalEphemeris" ).c_str( ) )
            .def( "body_fixed_to_inertial_rotation",
                  &te::RotationalEphemeris::getRotationMatrixToBaseFrame,
                  py::arg( "time" ),
                  get_docstring( "RotationalEphemeris.body_fixed_to_inertial_rotation" ).c_str( ) )
            .def( "time_derivative_body_fixed_to_inertial_rotation",
                  &te::RotationalEphemeris::getDerivativeOfRotationToBaseFrame,
                  py::arg( "time" ),
                  get_docstring( "RotationalEphemeris.time_derivative_body_fixed_to_inertial_rotation" ).c_str( ) )
            .def( "inertial_to_body_fixed_rotation",
                  &te::RotationalEphemeris::getRotationMatrixToTargetFrame,
                  py::arg( "time" ),
                  get_docstring( "RotationalEphemeris.inertial_to_body_fixed_rotation" ).c_str( ) )
            .def( "time_derivative_inertial_to_body_fixed_rotation",
                  &te::RotationalEphemeris::getDerivativeOfRotationToTargetFrame,
                  py::arg( "time" ),
                  get_docstring( "RotationalEphemeris.time_derivative_inertial_to_body_fixed_rotation" ).c_str( ) )
            .def( "angular_velocity_in_body_fixed_frame",
                  &te::RotationalEphemeris::getRotationalVelocityVectorInTargetFrame,
                  py::arg( "time" ),
                  get_docstring( "RotationalEphemeris.angular_velocity_in_body_fixed_frame" ).c_str( ) )
            .def( "angular_velocity_in_inertial_frame",
                  &te::RotationalEphemeris::getRotationalVelocityVectorInBaseFrame,
                  py::arg( "time" ),
                  get_docstring( "RotationalEphemeris.angular_velocity_in_inertial_frame" ).c_str( ) )
            .def_property_readonly( "body_fixed_frame_name",
                                    &te::RotationalEphemeris::getTargetFrameOrientation,
                                    get_docstring( "RotationalEphemeris.body_fixed_frame_name" ).c_str( ) )
            .def_property_readonly( "inertial_frame_name",
                                    &te::RotationalEphemeris::getBaseFrameOrientation,
                                    get_docstring( "RotationalEphemeris.inertial_frame_name" ).c_str( ) );

    m.def( "transform_to_inertial_orientation",
           &te::transformStateToInertialOrientation< double, double >,
           py::arg( "state_in_body_fixed_frame" ),
           py::arg( "current_time" ),
           py::arg( "rotational_ephemeris" ) );

    py::class_< te::LongitudeLibrationCalculator, std::shared_ptr< te::LongitudeLibrationCalculator > >( m,
                                                                                                         "LongitudeLibrationCalculator" );

    py::class_< te::DirectLongitudeLibrationCalculator,
                std::shared_ptr< te::DirectLongitudeLibrationCalculator >,
                te::LongitudeLibrationCalculator >( m, "DirectLongitudeLibrationCalculator" )
            .def( py::init< const double >( ), py::arg( "scaled_libration_amplitude" ) );

    py::class_< te::SynchronousRotationalEphemeris, std::shared_ptr< te::SynchronousRotationalEphemeris >, te::RotationalEphemeris >(
            m, "SynchronousRotationalEphemeris" )
            .def_property( "libration_calculator",
                           &te::SynchronousRotationalEphemeris::getLongitudeLibrationCalculator,
                           &te::SynchronousRotationalEphemeris::setLibrationCalculation );

    py::class_< te::AerodynamicAngleRotationalEphemeris,
                std::shared_ptr< te::AerodynamicAngleRotationalEphemeris >,
                te::RotationalEphemeris >( m, "AerodynamicAngleRotationalEphemeris" )
            .def( "reset_aerodynamic_angle_function", &te::AerodynamicAngleRotationalEphemeris::setAerodynamicAngleFunction );

    py::class_< te::GcrsToItrsRotationModel, std::shared_ptr< te::GcrsToItrsRotationModel >, te::RotationalEphemeris >(
            m, "GcrsToItrsRotationModel" );

    py::class_< te::DirectionBasedRotationalEphemeris, std::shared_ptr< te::DirectionBasedRotationalEphemeris >, te::RotationalEphemeris >(
            m, "CustomInertialDirectionBasedRotationalEphemeris" )
            .def_property_readonly( "inertial_body_axis_calculator",
                                    &te::DirectionBasedRotationalEphemeris::getInertialBodyAxisDirectionCalculator );

    py::class_< te::InertialBodyFixedDirectionCalculator, std::shared_ptr< te::InertialBodyFixedDirectionCalculator > >(
            m, "InertialBodyFixedDirectionCalculator" );

    py::class_< te::CustomBodyFixedDirectionCalculator,
                std::shared_ptr< te::CustomBodyFixedDirectionCalculator >,
                te::InertialBodyFixedDirectionCalculator >( m, "CustomBodyFixedDirectionCalculator" )
            .def_property( "inertial_body_axis_direction_function",
                           &te::CustomBodyFixedDirectionCalculator::getInertialBodyAxisDirectionFunction,
                           &te::CustomBodyFixedDirectionCalculator::resetInertialBodyAxisDirectionFunction );

    /*!
     **************   GRAVITY FIELD  ******************
     */

    py::class_< tg::GravityFieldModel, std::shared_ptr< tg::GravityFieldModel > >(
            m, "GravityFieldModel", get_docstring( "AtmosphereModel" ).c_str( ) )
            .def( py::init< const double, const std::function< void( ) > >( ),
                  py::arg( "gravitational_parameter" ),
                  py::arg( "update_inertia_tensor" ) = std::function< void( ) >( )  // <pybind11/functional.h>
                  )
            .def( "get_gravitational_parameter", &tg::GravityFieldModel::getGravitationalParameter )
            .def_property( "gravitational_parameter",
                           &tg::GravityFieldModel::getGravitationalParameter,
                           &tg::GravityFieldModel::resetGravitationalParameter );

    py::class_< tg::SphericalHarmonicsGravityField, std::shared_ptr< tg::SphericalHarmonicsGravityField >, tg::GravityFieldModel >(
            m, "SphericalHarmonicsGravityField" )
            .def_property_readonly( "reference_radius", &tg::SphericalHarmonicsGravityField::getReferenceRadius )
            .def_property_readonly( "maximum_degree", &tg::SphericalHarmonicsGravityField::getDegreeOfExpansion )
            .def_property_readonly( "maximum_order", &tg::SphericalHarmonicsGravityField::getOrderOfExpansion )
            .def_property( "cosine_coefficients",
                           &tg::SphericalHarmonicsGravityField::getCosineCoefficients,
                           &tg::SphericalHarmonicsGravityField::setCosineCoefficients )
            .def_property( "sine_coefficients",
                           &tg::SphericalHarmonicsGravityField::getSineCoefficients,
                           &tg::SphericalHarmonicsGravityField::setSineCoefficients );

    py::class_< tg::PolyhedronGravityField, std::shared_ptr< tg::PolyhedronGravityField >, tg::GravityFieldModel >(
            m, "PolyhedronGravityField" )
            .def_property_readonly( "volume", &tg::PolyhedronGravityField::getVolume )
            .def_property_readonly( "vertices_coordinates", &tg::PolyhedronGravityField::getVerticesCoordinates )
            .def_property_readonly( "vertices_defining_each_facet", &tg::PolyhedronGravityField::getVerticesDefiningEachFacet );

    /*!
     **************   SHAPE MODELS  ******************
     */

    py::class_< tba::BodyShapeModel, std::shared_ptr< tba::BodyShapeModel > >( m, "ShapeModel", get_docstring( "ShapeModel" ).c_str( ) )
            .def( "get_average_radius", &tba::BodyShapeModel::getAverageRadius )
            .def_property_readonly( "average_radius", &tba::BodyShapeModel::getAverageRadius );

    /*!
     **************   GROUND STATION FUNCTIONALITY  ******************
     */

    py::class_< tgs::GroundStationState, std::shared_ptr< tgs::GroundStationState > >( m, "GroundStationState" )
            .def( "get_cartesian_state",
                  &tgs::GroundStationState::getCartesianStateInTime,
                  py::arg( "seconds_since_epoch" ),
                  py::arg( "target_frame_origin" ) )
            .def( "get_cartesian_position",
                  &tgs::GroundStationState::getCartesianPositionInTime,
                  py::arg( "seconds_since_epoch" ),
                  py::arg( "target_frame_origin" ) )
            .def_property_readonly( "cartesian_positon_at_reference_epoch", &tgs::GroundStationState::getNominalCartesianPosition )
            .def_property_readonly( "spherical_positon_at_reference_epoch", &tgs::GroundStationState::getNominalSphericalPosition )
            .def_property_readonly( "geodetic_positon_at_reference_epoch", &tgs::GroundStationState::getNominalGeodeticPosition )
            .def_property_readonly( "rotation_matrix_body_fixed_to_topocentric",
                                    &tgs::GroundStationState::getRotationMatrixFromBodyFixedToTopocentricFrame );

    py::class_< tgs::GroundStation, std::shared_ptr< tgs::GroundStation > >( m, "GroundStation" )
            .def( "set_transmitting_frequency_calculator",
                  &tgs::GroundStation::setTransmittingFrequencyCalculator,
                  py::arg( "transmitting_frequency_calculator" ) )
            .def( "set_water_vapor_partial_pressure_function",
                  &tgs::GroundStation::setWaterVaporPartialPressureFunction,
                  py::arg( "water_vapor_partial_pressure_function" ) )
            .def( "set_temperature_function", &tgs::GroundStation::setTemperatureFunction, py::arg( "temperature_function" ) )
            .def( "set_pressure_function", &tgs::GroundStation::setPressureFunction, py::arg( "pressure_function" ) )
            .def( "set_relative_humidity_function",
                  &tgs::GroundStation::setRelativeHumidityFunction,
                  py::arg( "relative_humidity_function" ) )
            .def_property_readonly( "temperature_function", &tgs::GroundStation::getTemperatureFunction )
            .def_property_readonly( "pressure_function", &tgs::GroundStation::getPressureFunction )
            .def_property_readonly( "relative_humidity_function", &tgs::GroundStation::getRelativeHumidityFunction )
            .def_property_readonly( "pointing_angles_calculator", &tgs::GroundStation::getPointingAnglesCalculator )
            .def_property_readonly( "station_state", &tgs::GroundStation::getNominalStationState )
            .def( "set_timing_system", &tgs::GroundStation::setTimingSystem, py::arg( "timing_system" ) );

    py::class_< tgs::StationFrequencyInterpolator, std::shared_ptr< tgs::StationFrequencyInterpolator > >(
            m, "StationFrequencyInterpolator", get_docstring( "StationFrequencyInterpolator" ).c_str( ) );

    py::class_< tgs::ConstantFrequencyInterpolator,
                std::shared_ptr< tgs::ConstantFrequencyInterpolator >,
                tgs::StationFrequencyInterpolator >( m, "ConstantFrequencyInterpolator" )
            .def( py::init< double >( ), py::arg( "frequency" ) );

    py::class_< tgs::PointingAnglesCalculator, std::shared_ptr< tgs::PointingAnglesCalculator > >( m, "PointingAnglesCalculator" )
            .def( "calculate_elevation_angle",
                  py::overload_cast< const Eigen::Vector3d &, const double >(
                          &tgs::PointingAnglesCalculator::calculateElevationAngleFromInertialVector ),
                  py::arg( "inertial_vector_to_target" ),
                  py::arg( "time" ) )
            .def( "calculate_azimuth_angle",
                  py::overload_cast< const Eigen::Vector3d &, const double >(
                          &tgs::PointingAnglesCalculator::calculateAzimuthAngleFromInertialVector ),
                  py::arg( "inertial_vector_to_target" ),
                  py::arg( "time" ) )
            .def( "convert_inertial_vector_to_topocentric",
                  &tgs::PointingAnglesCalculator::convertVectorFromInertialToTopocentricFrame,
                  py::arg( "inertial_vector" ),
                  py::arg( "time" ) );

    /*!
     **************   BODY OBJECTS AND ASSOCIATED FUNCTIONALITY  ******************
     */

    py::class_< tss::Body, std::shared_ptr< tss::Body > >( m, "Body", get_docstring( "Body" ).c_str( ) )
            .def_property(
                    "ephemeris_frame_to_base_frame", &tss::Body::getEphemerisFrameToBaseFrame, &tss::Body::setEphemerisFrameToBaseFrame )
            .def_property_readonly( "state", &tss::Body::getState, get_docstring( "Body.state" ).c_str( ) )
            .def_property_readonly( "position", &tss::Body::getPosition, get_docstring( "Body.position" ).c_str( ) )
            .def_property_readonly( "velocity", &tss::Body::getVelocity, get_docstring( "Body.velocity" ).c_str( ) )
            .def_property_readonly( "inertial_to_body_fixed_frame",
                                    &tss::Body::getCurrentRotationMatrixToLocalFrame,
                                    get_docstring( "Body.inertial_to_body_fixed_frame" ).c_str( ) )
            .def_property_readonly( "body_fixed_to_inertial_frame",
                                    &tss::Body::getCurrentRotationMatrixToGlobalFrame,
                                    get_docstring( "Body.body_fixed_to_inertial_frame" ).c_str( ) )
            .def_property_readonly( "inertial_to_body_fixed_frame_derivative",
                                    &tss::Body::getCurrentRotationMatrixDerivativeToLocalFrame,
                                    get_docstring( "Body.inertial_to_body_fixed_frame_derivative" ).c_str( ) )
            .def_property_readonly( "body_fixed_to_inertial_frame_derivative",
                                    &tss::Body::getCurrentRotationMatrixDerivativeToGlobalFrame,
                                    get_docstring( "Body.body_fixed_to_inertial_frame_derivative" ).c_str( ) )
            .def_property_readonly( "inertial_angular_velocity",
                                    &tss::Body::getCurrentAngularVelocityVectorInGlobalFrame,
                                    get_docstring( "Body.inertial_angular_velocity" ).c_str( ) )
            .def_property_readonly( "body_fixed_angular_velocity",
                                    &tss::Body::getCurrentAngularVelocityVectorInLocalFrame,
                                    get_docstring( "Body.body_fixed_angular_velocity" ).c_str( ) )
            .def_property( "mass", &tss::Body::getBodyMass, &tss::Body::setConstantBodyMass, get_docstring( "Body.mass" ).c_str( ) )
            .def( "set_constant_mass", &tss::Body::setConstantBodyMass, py::arg( "mass" ) )
            .def_property( "inertia_tensor",
                           &tss::Body::getBodyInertiaTensor,
                           py::overload_cast< const Eigen::Matrix3d & >( &tss::Body::setBodyInertiaTensor ),
                           get_docstring( "Body.inertia_tensor" ).c_str( ) )
            .def( "state_in_base_frame_from_ephemeris",
                  &tss::Body::getStateInBaseFrameFromEphemeris< STATE_SCALAR_TYPE, TIME_TYPE >,
                  py::arg( "time" ) )
            .def_property( "ephemeris", &tss::Body::getEphemeris, &tss::Body::setEphemeris, get_docstring( "Body.ephemeris" ).c_str( ) )
            .def_property( "atmosphere_model",
                           &tss::Body::getAtmosphereModel,
                           &tss::Body::setAtmosphereModel,
                           get_docstring( "Body.atmosphere_model" ).c_str( ) )
            .def_property(
                    "shape_model", &tss::Body::getShapeModel, &tss::Body::setShapeModel, get_docstring( "Body.shape_model" ).c_str( ) )
            .def_property( "gravity_field_model",
                           &tss::Body::getGravityFieldModel,
                           &tss::Body::setGravityFieldModel,
                           get_docstring( "Body.gravity_field_model" ).c_str( ) )
            .def_property( "aerodynamic_coefficient_interface",
                           &tss::Body::getAerodynamicCoefficientInterface,
                           &tss::Body::setAerodynamicCoefficientInterface,
                           get_docstring( "Body.aerodynamic_coefficient_interface" ).c_str( ) )
            .def_property( "flight_conditions",
                           &tss::Body::getFlightConditions,
                           &tss::Body::setFlightConditions,
                           get_docstring( "Body.flight_conditions" ).c_str( ) )
            .def_property( "rotation_model",
                           &tss::Body::getRotationalEphemeris,
                           &tss::Body::setRotationalEphemeris,
                           get_docstring( "Body.rotation_model" ).c_str( ) )
            .def_property( "system_models",
                           &tss::Body::getVehicleSystems,
                           &tss::Body::setVehicleSystems,
                           get_docstring( "Body.system_models" ).c_str( ) )
            .def_property( "rigid_body_properties",
                           &tss::Body::getMassProperties,
                           &tss::Body::setMassProperties,
                           get_docstring( "Body.rigid_body_properties" ).c_str( ) )
            .def_property_readonly( "gravitational_parameter",
                                    &tss::Body::getGravitationalParameter,
                                    get_docstring( "Body.gravitational_parameter" ).c_str( ) )
            .def( "get_ground_station",
                  &tss::Body::getGroundStation,
                  py::arg( "station_name" ),
                  get_docstring( "Body.get_ground_station" ).c_str( ) )
            .def_property_readonly(
                    "ground_station_list", &tss::Body::getGroundStationMap, get_docstring( "Body.ground_station_list" ).c_str( ) );

    py::class_< tss::SystemOfBodies, std::shared_ptr< tss::SystemOfBodies > >(
            m, "SystemOfBodies", get_docstring( "SystemOfBodies" ).c_str( ) )
            .def( "get", &tss::SystemOfBodies::getBody, py::arg( "body_name" ), get_docstring( "SystemOfBodies.get" ).c_str( ) )
            .def( "get_body", &tss::SystemOfBodies::getBody, py::arg( "body_name" ), get_docstring( "SystemOfBodies.get_body" ).c_str( ) )
            .def( "create_empty_body",
                  &tss::SystemOfBodies::createEmptyBody< STATE_SCALAR_TYPE, TIME_TYPE >,
                  py::arg( "body_name" ),
                  py::arg( "process_body" ) = 1,
                  get_docstring( "SystemOfBodies.create_empty_body" ).c_str( ) )
            .def( "does_body_exist",
                  &tss::SystemOfBodies::doesBodyExist,
                  py::arg( "body_name" ),
                  get_docstring( "SystemOfBodies.does_body_exist" ).c_str( ) )
            .def( "list_of_bodies", &tss::SystemOfBodies::getListOfBodies, get_docstring( "SystemOfBodies.list_of_bodies" ).c_str( ) )
            //            .def("get_body_dict", &tss::SystemOfBodies::getMap,
            //                 get_docstring("SystemOfBodies.get_body_dict").c_str())
            .def( "add_body",
                  &tss::SystemOfBodies::addBody< STATE_SCALAR_TYPE, TIME_TYPE >,
                  py::arg( "body_to_add" ),
                  py::arg( "body_name" ),
                  py::arg( "process_body" ) = 1,
                  get_docstring( "SystemOfBodies.add_body" ).c_str( ) )
            .def( "remove_body",
                  &tss::SystemOfBodies::deleteBody,
                  py::arg( "body_name" ),
                  get_docstring( "SystemOfBodies.remove_body" ).c_str( ) )
            .def( "global_frame_orientation",
                  &tss::SystemOfBodies::getFrameOrientation,
                  get_docstring( "SystemOfBodies.global_frame_orientation" ).c_str( ) )
            .def( "global_frame_origin",
                  &tss::SystemOfBodies::getFrameOrigin,
                  get_docstring( "SystemOfBodies.global_frame_origin" ).c_str( ) );

    //            .def_property_readonly("number_of_bodies", &tss::SystemOfBodies::getNumberOfBodies,
    //                                   get_docstring("number_of_bodies").c_str() );

    /*!
     **************   SUPPORTING FUNCTIONS USED ENVIRONMENT MODELS  ******************
     */
}
}  // namespace environment
}  // namespace numerical_simulation
}  // namespace tudatpy

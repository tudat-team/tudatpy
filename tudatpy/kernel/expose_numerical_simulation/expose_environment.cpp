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
#include <tudat/basics/deprecationWarnings.h>

#include "tudatpy/docstrings.h"

#include <tudat/astro/aerodynamics.h>
#include <tudat/astro/ephemerides.h>
#include <tudat/astro/gravitation.h>
#include "tudat/astro/ground_stations/groundStation.h"
#include "tudat/simulation/environment_setup/body.h"


#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

//namespace py = pybind11;
//namespace tba = tudat::basic_astrodynamics;
//namespace tss = tudat::simulation_setup;
//namespace tp = tudat::propagators;
//namespace tinterp = tudat::interpolators;
//namespace te = tudat::ephemerides;
//namespace tni = tudat::numerical_integrators;
//namespace trf = tudat::reference_frames;
//namespace tmrf = tudat::root_finders;


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
    std::vector< boost::multi_array< Eigen::Vector3d, 2 > > meshPoints =
            localInclinationAnalysis->getMeshPoints( );
    std::vector< boost::multi_array< Eigen::Vector3d, 2 > > meshSurfaceNormals =
            localInclinationAnalysis->getPanelSurfaceNormals( );


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
                //                pressureCoefficientsList[ counter ] = ( Eigen::Vector1d( ) << pressureCoefficients[ i ][ j ][ k ] ).finished( );
                counter++;
            }
        }
    }

    return std::make_pair( meshPointsList, surfaceNormalsList );
}

}

}

namespace tudatpy {
namespace numerical_simulation {
namespace environment {

void expose_environment(py::module &m) {

    py::enum_<ta::AerodynamicCoefficientsIndependentVariables>(m, "AerodynamicCoefficientsIndependentVariables",
                                                               get_docstring("AerodynamicCoefficientsIndependentVariables").c_str())
            .value("mach_number_dependent", ta::AerodynamicCoefficientsIndependentVariables::mach_number_dependent,
                        get_docstring("AerodynamicCoefficientsIndependentVariables.mach_number_dependent").c_str())
            .value("angle_of_attack_dependent", ta::AerodynamicCoefficientsIndependentVariables::angle_of_attack_dependent,
                        get_docstring("AerodynamicCoefficientsIndependentVariables.angle_of_attack_dependent").c_str())
            .value("sideslip_angle_dependent", ta::AerodynamicCoefficientsIndependentVariables::angle_of_sideslip_dependent,
                        get_docstring("AerodynamicCoefficientsIndependentVariables.sideslip_angle_dependent").c_str())
            .value("altitude_dependent", ta::AerodynamicCoefficientsIndependentVariables::altitude_dependent,
                        get_docstring("AerodynamicCoefficientsIndependentVariables.altitude_dependent").c_str())
            .value("time_dependent", ta::AerodynamicCoefficientsIndependentVariables::time_dependent,
                        get_docstring("AerodynamicCoefficientsIndependentVariables.time_dependent").c_str())
            .value("control_surface_deflection_dependent", ta::AerodynamicCoefficientsIndependentVariables::control_surface_deflection_dependent,
                        get_docstring("AerodynamicCoefficientsIndependentVariables.control_surface_deflection_dependent").c_str())
            .value("undefined_independent_variable", ta::AerodynamicCoefficientsIndependentVariables::undefined_independent_variable,
                        get_docstring("AerodynamicCoefficientsIndependentVariables.undefined_independent_variable").c_str())
            .export_values();

    py::class_<ta::AerodynamicCoefficientInterface,
            std::shared_ptr<ta::AerodynamicCoefficientInterface>>(m, "AerodynamicCoefficientInterface" )
            .def_property_readonly("reference_area", &ta::AerodynamicCoefficientInterface::getReferenceArea )
            .def_property_readonly("current_force_coefficients", &ta::AerodynamicCoefficientInterface::getCurrentForceCoefficients )
            .def_property_readonly("current_moment_coefficients", &ta::AerodynamicCoefficientInterface::getCurrentMomentCoefficients )
            .def_property_readonly("current_coefficients", &ta::AerodynamicCoefficientInterface::getCurrentAerodynamicCoefficients )
            .def_property_readonly("independent_variable_names", &ta::AerodynamicCoefficientInterface::getIndependentVariableNames )
            .def_property_readonly("control_surface_independent_variable_names", &ta::AerodynamicCoefficientInterface::getControlSurfaceIndependentVariables )
            .def_property_readonly("current_control_surface_free_force_coefficients", &ta::AerodynamicCoefficientInterface::getCurrentControlSurfaceFreeForceCoefficients )
            .def_property_readonly("current_control_surface_free_moment_coefficients", &ta::AerodynamicCoefficientInterface::getCurrentControlSurfaceFreeMomentCoefficients )
            .def("current_control_surface_force_coefficient_increment", &ta::AerodynamicCoefficientInterface::getCurrentForceCoefficientIncrement,
                 py::arg( "control_surface_name" ) )
            .def("current_control_surface_moment_coefficient_increment", &ta::AerodynamicCoefficientInterface::getCurrentMomentCoefficientIncrement,
                 py::arg( "control_surface_name" ) )
            .def("set_control_surface_increments", &ta::AerodynamicCoefficientInterface::setControlSurfaceIncrements,
                 py::arg( "control_surface_list" ) )
            .def("update_coefficients", &ta::AerodynamicCoefficientInterface::updateCurrentCoefficients,
                 py::arg( "independent_variables" ),
                 py::arg( "time") )
            .def("update_full_coefficients", &ta::AerodynamicCoefficientInterface::updateFullCurrentCoefficients,
                 py::arg( "independent_variables" ),
                 py::arg( "control_surface_independent_variables" ),
                 py::arg( "time") );

    py::class_<ta::ControlSurfaceIncrementAerodynamicInterface,
            std::shared_ptr<ta::ControlSurfaceIncrementAerodynamicInterface>>(
                m, "ControlSurfaceIncrementAerodynamicInterface", "<no_doc, only_dec>");

    py::class_<ta::CustomControlSurfaceIncrementAerodynamicInterface,
            std::shared_ptr<ta::CustomControlSurfaceIncrementAerodynamicInterface>,
            ta::ControlSurfaceIncrementAerodynamicInterface>(
                m, "CustomControlSurfaceIncrementAerodynamicInterface", "<no_doc, only_dec>")
            .def(py::init<
                 const std::function< Eigen::Vector6d( const std::vector< double >& ) >,
                 const std::vector< ta::AerodynamicCoefficientsIndependentVariables > >(),
                 py::arg("coefficient_function"),
                 py::arg("independent_variable_names") );

    py::class_<ta::AerodynamicCoefficientGenerator<3, 6>,
            std::shared_ptr<ta::AerodynamicCoefficientGenerator<3, 6>>,
            ta::AerodynamicCoefficientInterface>(m, "AerodynamicCoefficientGenerator36", "<no_doc, only_dec>");


    m.def("get_default_local_inclination_mach_points", &ta::getDefaultHypersonicLocalInclinationMachPoints,
          py::arg( "mach_regime" ) = "Full" );

    m.def("get_default_local_inclination_angle_of_attack_points", &ta::getDefaultHypersonicLocalInclinationAngleOfAttackPoints );

    m.def("get_default_local_inclination_sideslip_angle_points", &ta::getDefaultHypersonicLocalInclinationAngleOfSideslipPoints );

    py::class_<ta::HypersonicLocalInclinationAnalysis,
            std::shared_ptr<ta::HypersonicLocalInclinationAnalysis>,
            ta::AerodynamicCoefficientGenerator<3, 6>>(m, "HypersonicLocalInclinationAnalysis" )
            .def(py::init<
                 const std::vector< std::vector< double > >&,
                 const std::shared_ptr< tudat::SurfaceGeometry >,
                 const std::vector< int >&,
                 const std::vector< int >&,
                 const std::vector< bool >&,
                 const std::vector< std::vector< int > >&,
                 const double,
                 const double,
                 const Eigen::Vector3d&,
                 const bool >(),
                 py::arg("independent_variable_points"),
                 py::arg("body_shape"),
                 py::arg("number_of_lines"),
                 py::arg("number_of_points"),
                 py::arg("invert_orders"),
                 py::arg("selected_methods"),
                 py::arg("reference_area"),
                 py::arg("reference_length"),
                 py::arg("moment_reference_point"),
                 py::arg("save_pressure_coefficients") = false,
                 get_docstring("HypersonicLocalInclinationAnalysis.ctor").c_str())
            .def("clear_data",
                 &ta::HypersonicLocalInclinationAnalysis::clearData );


    m.def("save_vehicle_mesh_to_file", &ta::saveVehicleMeshToFile,
          py::arg( "local_inclination_analysis_object" ),
          py::arg( "output_directory" ),
          py::arg( "output_file_prefix" ) = "",
          get_docstring("save_vehicle_mesh_to_file").c_str() );

    m.def("get_local_inclination_total_vehicle_area", &ta::getTotalSurfaceArea,
          py::arg( "local_inclination_analysis_object" ) );

    m.def("get_local_inclination_mesh", &ta::getVehicleMesh,
          py::arg( "local_inclination_analysis_object" ) );




    py::class_<tsm::VehicleSystems,
            std::shared_ptr<tsm::VehicleSystems>>(m, "VehicleSystems" )
            .def(py::init< >() )
            .def("set_control_surface_deflection",
                 &tsm::VehicleSystems::setCurrentControlSurfaceDeflection,
                 py::arg("control_surface_id"),
                 py::arg("deflection_angle"));

    /*!
     **************   FLIGHT CONDITIONS AND ASSOCIATED FUNCTIONALITY  ******************
     */

    py::enum_<trf::AerodynamicsReferenceFrameAngles>(m, "AerodynamicsReferenceFrameAngles")
            .value("latitude_angle", trf::AerodynamicsReferenceFrameAngles::latitude_angle)
            .value("longitude_angle", trf::AerodynamicsReferenceFrameAngles::longitude_angle)
            .value("heading_angle", trf::AerodynamicsReferenceFrameAngles::heading_angle)
            .value("flight_path_angle", trf::AerodynamicsReferenceFrameAngles::flight_path_angle)
            .value("angle_of_attack", trf::AerodynamicsReferenceFrameAngles::angle_of_attack)
            .value("angle_of_sideslip", trf::AerodynamicsReferenceFrameAngles::angle_of_sideslip)
            .value("bank_angle", trf::AerodynamicsReferenceFrameAngles::bank_angle)
            .export_values();

    py::enum_<trf::AerodynamicsReferenceFrames>(m, "AerodynamicsReferenceFrames",
                                                get_docstring("AerodynamicsReferenceFrames").c_str())
            .value("inertial_frame", trf::AerodynamicsReferenceFrames::inertial_frame,
                   get_docstring("AerodynamicsReferenceFrames.inertial_frame").c_str())
            .value("corotating_frame", trf::AerodynamicsReferenceFrames::corotating_frame,
                   get_docstring("AerodynamicsReferenceFrames.corotating_frame").c_str())
            .value("vertical_frame", trf::AerodynamicsReferenceFrames::vertical_frame,
                   get_docstring("AerodynamicsReferenceFrames.vertical_frame").c_str())
            .value("trajectory_frame", trf::AerodynamicsReferenceFrames::trajectory_frame,
                   get_docstring("AerodynamicsReferenceFrames.trajectory_frame").c_str())
            .value("aerodynamic_frame", trf::AerodynamicsReferenceFrames::aerodynamic_frame,
                   get_docstring("AerodynamicsReferenceFrames.aerodynamic_frame").c_str())
            .value("body_frame", trf::AerodynamicsReferenceFrames::body_frame,
                   get_docstring("AerodynamicsReferenceFrames.body_frame").c_str())
            .export_values();

    py::class_<trf::AerodynamicAngleCalculator,
            std::shared_ptr<trf::AerodynamicAngleCalculator>>(m, "AerodynamicAngleCalculator", get_docstring("AerodynamicAngleCalculator").c_str())
            .def("get_rotation_matrix_between_frames",
                 &trf::AerodynamicAngleCalculator::getRotationMatrixBetweenFrames,
                 py::arg("original_frame"),
                 py::arg("target_frame"),
                 get_docstring("AerodynamicAngleCalculator.get_rotation_matrix_between_frames").c_str())
            .def("get_angle",
                 &trf::AerodynamicAngleCalculator::getAerodynamicAngle,
                 py::arg("angle_type"),
                 get_docstring("AerodynamicAngleCalculator.get_angle").c_str())
            // Function removed; error is shown
            .def("set_body_orientation_angles",
                 &trf::AerodynamicAngleCalculator::setOrientationAngleFunctionsRemoved2,
                 py::arg("angle_of_attack") = TUDAT_NAN,
                 py::arg("angle_of_sideslip") = TUDAT_NAN,
                 py::arg("bank_angle") = TUDAT_NAN,
                 py::arg("silence_warnings")=false)
            // Function removed; error is shown
            .def("set_body_orientation_angle_functions",
                 &trf::AerodynamicAngleCalculator::setOrientationAngleFunctionsRemoved1,
                 py::arg("angle_of_attack_function") = std::function<double()>(),       // <pybind11/functional.h>
                 py::arg("angle_of_sideslip_function") = std::function<double()>(),     // <pybind11/functional.h>
                 py::arg("bank_angle_function") = std::function<double()>(),            // <pybind11/functional.h>
                 py::arg("angle_update_function") = std::function<void(
                const double)>(),
                 py::arg("silence_warnings")=false);



    py::class_<ta::FlightConditions,
            std::shared_ptr<ta::FlightConditions>>(m, "FlightConditions", get_docstring("FlightConditions").c_str())
//            .def(py::init<
//                 const std::shared_ptr<tudat::basic_astrodynamics::BodyShapeModel>,
//                 const std::shared_ptr<tudat::reference_frames::AerodynamicAngleCalculator>>(),
//                 py::arg("shape_model"),
//                 py::arg("aerodynamic_angle_calculator") = std::shared_ptr< tr::AerodynamicAngleCalculator>())
            .def("update_conditions", &ta::FlightConditions::updateConditions, py::arg("current_time") )
            .def_property_readonly("aerodynamic_angle_calculator", &ta::FlightConditions::getAerodynamicAngleCalculator, get_docstring("FlightConditions.aerodynamic_angle_calculator").c_str())
            .def_property_readonly("longitude", &ta::FlightConditions::getCurrentLongitude, get_docstring("FlightConditions.longitude").c_str())
            .def_property_readonly("geodetic_latitude", &ta::FlightConditions::getCurrentGeodeticLatitude, get_docstring("FlightConditions.latitude").c_str())
            .def_property_readonly("time", &ta::FlightConditions::getCurrentTime, get_docstring("FlightConditions.time").c_str())
            .def_property_readonly("body_centered_body_fixed_state", &ta::FlightConditions::getCurrentBodyCenteredBodyFixedState, get_docstring("FlightConditions.body_centered_body_fixed_state").c_str())
            .def_property_readonly("altitude", &ta::FlightConditions::getCurrentAltitude, get_docstring("FlightConditions.time").c_str());

    py::class_<ta::AtmosphericFlightConditions,
            std::shared_ptr<ta::AtmosphericFlightConditions>,
            ta::FlightConditions>(m, "AtmosphericFlightConditions", get_docstring("AtmosphericFlightConditions").c_str())
            .def_property_readonly("density", &ta::AtmosphericFlightConditions::getCurrentDensity, get_docstring("AtmosphericFlightConditions.density").c_str())
            .def_property_readonly("temperature", &ta::AtmosphericFlightConditions::getCurrentFreestreamTemperature, get_docstring("AtmosphericFlightConditions.temperature").c_str())
            .def_property_readonly("dynamic_pressure", &ta::AtmosphericFlightConditions::getCurrentDynamicPressure, get_docstring("AtmosphericFlightConditions.dynamic_pressure").c_str())
            .def_property_readonly("pressure", &ta::AtmosphericFlightConditions::getCurrentPressure, get_docstring("AtmosphericFlightConditions.pressure").c_str())
            .def_property_readonly("airspeed", &ta::AtmosphericFlightConditions::getCurrentAirspeed, get_docstring("AtmosphericFlightConditions.airspeed").c_str())
            .def_property_readonly("mach_number", &ta::AtmosphericFlightConditions::getCurrentMachNumber, get_docstring("AtmosphericFlightConditions.mach_number").c_str())
            .def_property_readonly("airspeed_velocity", &ta::AtmosphericFlightConditions::getCurrentAirspeedBasedVelocity, get_docstring("AtmosphericFlightConditions.airspeed_velocity").c_str())
            .def_property_readonly("speed_of_sound", &ta::AtmosphericFlightConditions::getCurrentSpeedOfSound, get_docstring("AtmosphericFlightConditions.speed_of_sound").c_str())
            .def_property_readonly("aero_coefficient_independent_variables",
                                   &ta::AtmosphericFlightConditions::getAerodynamicCoefficientIndependentVariables, get_docstring("AtmosphericFlightConditions.aero_coefficient_independent_variables").c_str())
            .def_property_readonly("control_surface_aero_coefficient_independent_variables",
                                   &ta::AtmosphericFlightConditions::getControlSurfaceAerodynamicCoefficientIndependentVariables, get_docstring("AtmosphericFlightConditions.control_surface_aero_coefficient_independent_variables").c_str())
            .def_property_readonly("aerodynamic_coefficient_interface", &ta::AtmosphericFlightConditions::getAerodynamicCoefficientInterface, get_docstring("AtmosphericFlightConditions.aerodynamic_coefficient_interface").c_str());



    /*!
     **************   EPHEMERIDES  ******************
     */

    py::class_<te::Ephemeris, std::shared_ptr<te::Ephemeris>>(m, "Ephemeris")
            .def("cartesian_state", &te::Ephemeris::getCartesianState, py::arg("seconds_since_epoch") = 0.0)
            .def("cartesian_position", &te::Ephemeris::getCartesianPosition, py::arg("seconds_since_epoch") = 0.0)
            .def("cartesian_velocity", &te::Ephemeris::getCartesianVelocity, py::arg("seconds_since_epoch") = 0.0);


    py::class_<te::ConstantEphemeris,
            std::shared_ptr<te::ConstantEphemeris>,
            te::Ephemeris>(
                m, "ConstantEphemeris")
            .def(py::init<
                 const std::function<Eigen::Vector6d()>,//<pybind11/functional.h>,<pybind11/eigen.h>
                 const std::string &,
                 const std::string &>(),
                 py::arg("constant_state_function"),
                 py::arg("reference_frame_origin") = "SSB",
                 py::arg("reference_frame_orientation") = "ECLIPJ2000")
            .def(py::init<
                 const Eigen::Vector6d,//<pybind11/eigen.h>
                 const std::string &,
                 const std::string &>(),
                 py::arg("constant_state"),
                 py::arg("reference_frame_origin") = "SSB",
                 py::arg("reference_frame_orientation") = "ECLIPJ2000")
            .def("update_constant_state", &te::ConstantEphemeris::updateConstantState,
                 py::arg("new_state"));


    py::class_<te::KeplerEphemeris,
            std::shared_ptr<te::KeplerEphemeris>,
            te::Ephemeris>(
                m, "KeplerEphemeris");


    py::class_<te::MultiArcEphemeris, std::shared_ptr<te::MultiArcEphemeris>, te::Ephemeris>(m, "MultiArcEphemeris")
            .def(py::init<
                         const std::map<double, std::shared_ptr<te::Ephemeris> > &,
                         const std::string &,
                         const std::string &>(),
                 py::arg("single_arc_ephemerides"),
                 py::arg("reference_frame_origin") = "SSB",
                 py::arg("reference_frame_orientation") = "ECLIPJ2000");


    py::class_<te::TabulatedCartesianEphemeris< double, double >,
            std::shared_ptr<te::TabulatedCartesianEphemeris< double, double > >,
            te::Ephemeris>(m, "TabulatedEphemeris")
            .def_property("interpolator",
                          &te::TabulatedCartesianEphemeris< double, double >::getDynamicVectorSizeInterpolator,
                          py::overload_cast<
                          const std::shared_ptr< ti::OneDimensionalInterpolator
                          < double, Eigen::VectorXd > > >(
                              &te::TabulatedCartesianEphemeris< double, double >::resetInterpolator ) );


    py::class_<te::Tle, std::shared_ptr<te::Tle>>(m, "Tle")
            .def(py::init<//ctor 1
                 const std::string &>(),
                 py::arg("lines"))
            .def(py::init<//ctor 2
                 const std::string &,
                 const std::string &>(),
                 py::arg("line_1"),
                 py::arg("line_2"))
            .def("get_epoch", &te::Tle::getEpoch)
            .def("get_b_star", &te::Tle::getBStar)
            .def("get_epoch", &te::Tle::getEpoch)
            .def("get_inclination", &te::Tle::getInclination)
            .def("get_right_ascension", &te::Tle::getRightAscension)
            .def("get_eccentricity", &te::Tle::getEccentricity)
            .def("get_arg_of_perigee", &te::Tle::getArgOfPerigee)
            .def("get_mean_anomaly", &te::Tle::getMeanAnomaly)
            .def("get_mean_motion", &te::Tle::getMeanMotion);

    py::class_<te::TleEphemeris,
            std::shared_ptr<te::TleEphemeris>,
            te::Ephemeris>(m, "TleEphemeris")
            .def(py::init<
                 const std::string &,
                 const std::string &,
                 const std::shared_ptr<te::Tle>,
                 const bool>(),
                 py::arg("frame_origin") = "Earth",
                 py::arg("frame_orientation") = "J2000",
                 py::arg("tle") = nullptr,
                 py::arg("use_sdp") = false);

    /*!
     **************   ROTATION MODELS  ******************
     */


    py::class_<te::RotationalEphemeris,
            std::shared_ptr<te::RotationalEphemeris>>(m, "RotationalEphemeris")
            .def("body_fixed_to_inertial_rotation", &te::RotationalEphemeris::getRotationMatrixToBaseFrame,
                 py::arg( "time" ) )
            .def("time_derivative_body_fixed_to_inertial_rotation", &te::RotationalEphemeris::getDerivativeOfRotationToBaseFrame,
                 py::arg( "time" ) )
            .def("inertial_to_body_fixed_rotation", &te::RotationalEphemeris::getRotationMatrixToTargetFrame,
                 py::arg( "time" ) )
            .def("time_derivative_inertial_to_body_fixed_rotation", &te::RotationalEphemeris::getDerivativeOfRotationToTargetFrame,
                 py::arg( "time" ) )
            .def("angular_velocity_in_body_fixed_frame", &te::RotationalEphemeris::getRotationalVelocityVectorInTargetFrame,
                 py::arg( "time" ) )
            .def("angular_velocity_in_inertial_frame", &te::RotationalEphemeris::getRotationalVelocityVectorInBaseFrame,
                 py::arg( "time" ) )
            .def_property_readonly("body_fixed_frame_name", &te::RotationalEphemeris::getTargetFrameOrientation );


    m.def("transform_to_inertial_orientation",
          &te::transformStateToInertialOrientation<double, double>,
          py::arg("state_in_body_fixed_frame"),
          py::arg("current_time"),
          py::arg("rotational_ephemeris"));



    py::class_<te::LongitudeLibrationCalculator,
            std::shared_ptr<te::LongitudeLibrationCalculator>>(
                m, "LongitudeLibrationCalculator");

    py::class_<te::DirectLongitudeLibrationCalculator,
            std::shared_ptr<te::DirectLongitudeLibrationCalculator>,
            te::LongitudeLibrationCalculator>(
                m, "DirectLongitudeLibrationCalculator")
            .def(py::init< const double >(),
                 py::arg("scaled_libration_amplitude"));


    py::class_<te::SynchronousRotationalEphemeris,
            std::shared_ptr<te::SynchronousRotationalEphemeris>,
            te::RotationalEphemeris>(
                m, "SynchronousRotationalEphemeris")
            .def_property("libration_calculator",
                          &te::SynchronousRotationalEphemeris::getLongitudeLibrationCalculator,
                          &te::SynchronousRotationalEphemeris::setLibrationCalculation);

    py::class_<te::AerodynamicAngleRotationalEphemeris,
            std::shared_ptr<te::AerodynamicAngleRotationalEphemeris>,
            te::RotationalEphemeris>(
            m, "AerodynamicAngleRotationalEphemeris")
            .def("reset_aerodynamic_angle_function",
                          &te::AerodynamicAngleRotationalEphemeris::setAerodynamicAngleFunction );

    py::class_<te::GcrsToItrsRotationModel,
            std::shared_ptr<te::GcrsToItrsRotationModel>,
            te::RotationalEphemeris>(
                m, "GcrsToItrsRotationModel");


    /*!
     **************   GRAVITY FIELD  ******************
     */

    py::class_<tg::GravityFieldModel,
            std::shared_ptr<tg::GravityFieldModel>>(m, "GravityFieldModel")
            .def(py::init<
                 const double,
                 const std::function<void()>>(),
                 py::arg("gravitational_parameter"),
                 py::arg("update_inertia_tensor") = std::function<void()>()// <pybind11/functional.h>
            )
            .def("get_gravitational_parameter", &tg::GravityFieldModel::getGravitationalParameter)
            .def_property("gravitational_parameter", &tg::GravityFieldModel::getGravitationalParameter,
                          &tg::GravityFieldModel::resetGravitationalParameter);

    py::class_<tg::SphericalHarmonicsGravityField,
            std::shared_ptr<tg::SphericalHarmonicsGravityField >,
            tg::GravityFieldModel>(m, "SphericalHarmonicsGravityField")
            .def_property_readonly("reference_radius", &tg::SphericalHarmonicsGravityField::getReferenceRadius )
            .def_property_readonly("maximum_degree", &tg::SphericalHarmonicsGravityField::getDegreeOfExpansion )
            .def_property_readonly("maximum_order", &tg::SphericalHarmonicsGravityField::getOrderOfExpansion )
            .def_property("cosine_coefficients", &tg::SphericalHarmonicsGravityField::getCosineCoefficients,
                          &tg::SphericalHarmonicsGravityField::setCosineCoefficients)
            .def_property("sine_coefficients", &tg::SphericalHarmonicsGravityField::getSineCoefficients,
                          &tg::SphericalHarmonicsGravityField::setSineCoefficients);

    py::class_<tg::PolyhedronGravityField,
            std::shared_ptr<tg::PolyhedronGravityField >,
            tg::GravityFieldModel>(m, "PolyhedronGravityField")
            .def_property_readonly("volume", &tg::PolyhedronGravityField::getVolume )
            .def_property_readonly("vertices_coordinates", &tg::PolyhedronGravityField::getVerticesCoordinates )
            .def_property_readonly("vertices_defining_each_facet", &tg::PolyhedronGravityField::getVerticesDefiningEachFacet );

    /*!
     **************   SHAPE MODELS  ******************
     */

    py::class_<tba::BodyShapeModel,
            std::shared_ptr<tba::BodyShapeModel>>(m, "ShapeModel")
            .def("get_average_radius", &tba::BodyShapeModel::getAverageRadius)
            .def_property_readonly("average_radius", &tba::BodyShapeModel::getAverageRadius);


    /*!
     **************   GROUND STATION FUNCTIONALITY  ******************
     */


    py::class_<tgs::GroundStationState,
            std::shared_ptr<tgs::GroundStationState>>(m, "GroundStationState")
            .def("get_cartesian_state", &tgs::GroundStationState::getCartesianStateInTime,
                 py::arg( "seconds_since_epoch" ) )
            .def("get_cartesian_position", &tgs::GroundStationState::getCartesianPositionInTime,
                 py::arg( "seconds_since_epoch" ) )
            .def_property_readonly("cartesian_positon_at_reference_epoch", &tgs::GroundStationState::getNominalCartesianPosition )
            .def_property_readonly("spherical_positon_at_reference_epoch", &tgs::GroundStationState::getNominalSphericalPosition )
            .def_property_readonly("geodetic_positon_at_reference_epoch", &tgs::GroundStationState::getNominalGeodeticPosition )
            .def_property_readonly("rotation_matrix_body_fixed_to_topocentric", &tgs::GroundStationState::getRotationMatrixFromBodyFixedToTopocentricFrame );

    py::class_<tgs::GroundStation,
            std::shared_ptr<tgs::GroundStation>>(m, "GroundStation")
            .def_property_readonly("pointing_angles_calculator", &tgs::GroundStation::getPointingAnglesCalculator )
            .def_property_readonly("station_state", &tgs::GroundStation::getNominalStationState );


    py::class_<tgs::PointingAnglesCalculator,
            std::shared_ptr<tgs::PointingAnglesCalculator>>(m, "PointingAnglesCalculator")
            .def("calculate_elevation_angle", &tgs::PointingAnglesCalculator::calculateElevationAngle,
                 py::arg( "inertial_vector_to_target" ),
                 py::arg( "time" ) )
            .def("calculate_azimuth_angle", &tgs::PointingAnglesCalculator::calculateAzimuthAngle,
                 py::arg( "inertial_vector_to_target" ),
                 py::arg( "time" ) )
            .def("convert_inertial_vector_to_topocentric",
                 &tgs::PointingAnglesCalculator::convertVectorFromInertialToTopocentricFrame,
                 py::arg( "inertial_vector" ),
                 py::arg( "time" ) );


    /*!
     **************   BODY OBJECTS AND ASSOCIATED FUNCTIONALITY  ******************
     */

    py::class_<tss::Body, std::shared_ptr<tss::Body>>(m, "Body", get_docstring("Body").c_str())
            .def_property("ephemeris_frame_to_base_frame", &tss::Body::getEphemerisFrameToBaseFrame,
                          &tss::Body::setEphemerisFrameToBaseFrame)
            .def_property_readonly("state", &tss::Body::getState, get_docstring("Body.state").c_str())
            .def_property_readonly("position", &tss::Body::getPosition, get_docstring("Body.position").c_str())
            .def_property_readonly("velocity", &tss::Body::getVelocity, get_docstring("Body.velocity").c_str())
            .def_property_readonly("inertial_to_body_fixed_frame", &tss::Body::getCurrentRotationMatrixToLocalFrame, get_docstring("Body.inertial_to_body_fixed_frame").c_str())
            .def_property_readonly("body_fixed_to_inertial_frame", &tss::Body::getCurrentRotationMatrixToGlobalFrame, get_docstring("Body.body_fixed_to_inertial_frame").c_str())
            .def_property_readonly("inertial_to_body_fixed_frame_derivative", &tss::Body::getCurrentRotationMatrixDerivativeToLocalFrame, get_docstring("Body.inertial_to_body_fixed_frame_derivative").c_str())
            .def_property_readonly("body_fixed_to_inertial_frame_derivative", &tss::Body::getCurrentRotationMatrixDerivativeToGlobalFrame, get_docstring("Body.body_fixed_to_inertial_frame_derivative").c_str())
            .def_property_readonly("inertial_angular_velocity", &tss::Body::getCurrentAngularVelocityVectorInGlobalFrame, get_docstring("Body.inertial_angular_velocity").c_str())
            .def_property_readonly("body_fixed_angular_velocity", &tss::Body::getCurrentAngularVelocityVectorInLocalFrame, get_docstring("Body.body_fixed_angular_velocity").c_str())
            .def_property("mass", &tss::Body::getBodyMass, &tss::Body::setConstantBodyMass, get_docstring("Body.mass").c_str())
            .def("set_constant_mass", &tss::Body::setConstantBodyMass, py::arg( "mass" ) )
            .def_property("inertia_tensor", &tss::Body::getBodyInertiaTensor,
                          py::overload_cast<const Eigen::Matrix3d &>(
                              &tss::Body::setBodyInertiaTensor))
            .def("state_in_base_frame_from_ephemeris",
                 &tss::Body::getStateInBaseFrameFromEphemeris<double, double>, py::arg("time"))
            .def_property_readonly("ephemeris", &tss::Body::getEphemeris, get_docstring("Body.ephemeris").c_str())
            .def_property("atmosphere_model", &tss::Body::getAtmosphereModel, &tss::Body::setAtmosphereModel, get_docstring("Body.atmosphere_model").c_str())
            .def_property("shape_model", &tss::Body::getShapeModel, &tss::Body::setShapeModel, get_docstring("Body.shape_model").c_str())
            .def_property("gravity_field_model", &tss::Body::getGravityFieldModel, &tss::Body::setGravityFieldModel, get_docstring("Body.gravity_field_model").c_str())
            .def_property("aerodynamic_coefficient_interface", &tss::Body::getAerodynamicCoefficientInterface,
                          &tss::Body::setAerodynamicCoefficientInterface, get_docstring("Body.aerodynamic_coefficient_interface").c_str())
            .def_property("flight_conditions", &tss::Body::getFlightConditions, &tss::Body::setFlightConditions)
            .def_property("rotation_model", &tss::Body::getRotationalEphemeris, &tss::Body::setRotationalEphemeris, get_docstring("Body.rotation_model").c_str())
            .def_property("system_models", &tss::Body::getVehicleSystems, &tss::Body::setVehicleSystems)
            .def_property_readonly("gravitational_parameter", &tss::Body::getGravitationalParameter, get_docstring("Body.gravitational_parameter").c_str())
            .def("get_ground_station", &tss::Body::getGroundStation, py::arg("station_name"))
            .def_property_readonly("ground_station_list", &tss::Body::getGroundStationMap );




    py::class_<tss::SystemOfBodies,
            std::shared_ptr<tss::SystemOfBodies> >(m, "SystemOfBodies", get_docstring("SystemOfBodies").c_str())
            .def("get", &tss::SystemOfBodies::getBody,
                 py::arg("body_name"),
                 get_docstring("SystemOfBodies.get").c_str())
            .def("get_body", &tss::SystemOfBodies::getBody,
                 py::arg("body_name"),
                 get_docstring("SystemOfBodies.get_body").c_str())
            .def("create_empty_body", &tss::SystemOfBodies::createEmptyBody,
                 py::arg("body_name"),
                 py::arg("process_body") = 1,
                 get_docstring("SystemOfBodies.create_empty_body").c_str())
            .def("add_body", &tss::SystemOfBodies::addBody,
                 py::arg("body_to_add"),
                 py::arg("body_name"),
                 py::arg("process_body") = 1,
                 get_docstring("SystemOfBodies.add_body").c_str())
            .def("remove_body", &tss::SystemOfBodies::deleteBody,
                 py::arg("body_name"),
                 get_docstring("SystemOfBodies.remove_body").c_str());
//            .def_property_readonly("number_of_bodies", &tss::SystemOfBodies::getNumberOfBodies,
//                                   get_docstring("number_of_bodies").c_str() );

    /*!
     **************   SUPPORTING FUNCTIONS USED ENVIRONMENT MODELS  ******************
     */


}
}// namespace environment
}// namespace numerical_simulation
}// namespace tudatpy

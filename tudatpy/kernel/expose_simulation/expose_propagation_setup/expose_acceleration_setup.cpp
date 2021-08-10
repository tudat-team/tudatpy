/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_acceleration_setup.h"

#include "tudatpy/docstrings.h"
#include <tudat/simulation/propagation_setup.h>

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tinterp = tudat::interpolators;
namespace te = tudat::ephemerides;
namespace tni = tudat::numerical_integrators;
namespace trf = tudat::reference_frames;
namespace tmrf = tudat::root_finders;

namespace tudatpy {
namespace simulation {
namespace propagation_setup {

void expose_acceleration_setup(py::module &m) {

    /*
     * This contains the addition of IntegratorSettings and AvailableIntegrators
     * and AvailableAccelerations which should be relocated in the tudat source.
     */

    py::enum_<tba::AvailableAcceleration>(m, "AvailableAcceleration")
            .value("undefined_acceleration_type", tba::AvailableAcceleration::undefined_acceleration)
            .value("point_mass_gravity_type", tba::AvailableAcceleration::point_mass_gravity)
            .value("central_gravity_type", tba::AvailableAcceleration::central_gravity)
            .value("aerodynamic_type", tba::AvailableAcceleration::aerodynamic)
            .value("cannonball_radiation_pressure_type", tba::AvailableAcceleration::cannon_ball_radiation_pressure)
            .value("spherical_harmonic_gravity_type", tba::AvailableAcceleration::spherical_harmonic_gravity)
            .value("mutual_spherical_harmonic_gravity_type", tba::AvailableAcceleration::mutual_spherical_harmonic_gravity)
            .value("third_body_point_mass_gravity_type", tba::AvailableAcceleration::third_body_point_mass_gravity)
            .value("third_body_central_gravity_type", tba::AvailableAcceleration::third_body_central_gravity)
            .value("third_body_spherical_harmonic_gravity_type", tba::AvailableAcceleration::third_body_spherical_harmonic_gravity)
            .value("third_body_mutual_spherical_harmonic_gravity_type", tba::AvailableAcceleration::third_body_mutual_spherical_harmonic_gravity)
            .value("thrust_acceleration_type", tba::AvailableAcceleration::thrust_acceleration)
            .value("relativistic_correction_acceleration_type", tba::AvailableAcceleration::relativistic_correction_acceleration)
            .value("empirical_acceleration_type", tba::AvailableAcceleration::empirical_acceleration)
            .value("direct_tidal_dissipation_in_central_body_acceleration_type", tba::AvailableAcceleration::direct_tidal_dissipation_in_central_body_acceleration)
            .value("direct_tidal_dissipation_in_orbiting_body_acceleration_type", tba::AvailableAcceleration::direct_tidal_dissipation_in_orbiting_body_acceleration)
            .value("panelled_radiation_pressure_acceleration_type", tba::AvailableAcceleration::panelled_radiation_pressure_acceleration)
            .value("momentum_wheel_desaturation_acceleration_type", tba::AvailableAcceleration::momentum_wheel_desaturation_acceleration)
            .value("solar_sail_acceleration_type", tba::AvailableAcceleration::solar_sail_acceleration)
            .export_values();

    py::enum_<tss::ThrustFrames>(m, "ThrustFrames")
            .value("unspecified_thrust_frame_type", tss::ThrustFrames::unspecified_thrust_frame)
            .value("inertial_thrust_frame_type", tss::ThrustFrames::inertial_thrust_frame)
            .value("lvlh_thrust_frame_type", tss::ThrustFrames::lvlh_thrust_frame)
            .export_values();

    //////////////////////////////////////////////////////////////////////////////
    // accelerationSettings.h
    //////////////////////////////////////////////////////////////////////////////

    // Unified interface functions for acceleration settings
    //  m.def("acceleration", &tss::acceleration, py::arg("acceleration_type"));
    m.def("point_mass_gravity", &tss::pointMassGravityAcceleration,
          get_docstring("point_mass_gravity").c_str());

    m.def("aerodynamic", &tss::aerodynamicAcceleration,
          get_docstring("aerodynamic").c_str());

    m.def("cannonball_radiation_pressure", &tss::cannonBallRadiationPressureAcceleration,
          get_docstring("cannonball_radiation_pressure").c_str());

    m.def("spherical_harmonic_gravity", &tss::sphericalHarmonicAcceleration,
          py::arg( "maximum_degree" ),
          py::arg( "maximum_order" ),
          get_docstring("spherical_harmonic_gravity").c_str());

    m.def("mutual_spherical_harmonic_gravity", &tss::mutualSphericalHarmonicAcceleration,
          py::arg( "maximum_degree_body_exerting" ),
          py::arg( "maximum_order_body_exerting" ),
          py::arg( "maximum_degree_body_undergoing" ),
          py::arg( "maximum_order_body_undergoing" ),
          py::arg( "maximum_degree_central_body" ) = 0,
          py::arg( "maximum_order_central_body" ) = 0,
          get_docstring("mutual_spherical_harmonic_gravity").c_str());

    // TODO: I replaced zeros with false, is it okay?
    m.def("relativistic_correction", &tss::relativisticAccelerationCorrection,
          py::arg( "use_schwarzschild" ) = false,
          py::arg( "use_lense_thirring" ) = false,
          py::arg( "use_de_sitter" ) = false,
          py::arg( "de_sitter_central_body" ) = "",
          py::arg( "lense_thirring_angular_momentum" ) = Eigen::Vector3d::Zero( ),
          get_docstring("relativistic_correction").c_str());

    m.def("empirical", &tss::empiricalAcceleration,
          py::arg( "constant_acceleration" ) = Eigen::Vector3d::Zero( ),
          py::arg( "sine_acceleration" ) = Eigen::Vector3d::Zero( ),
          py::arg( "cosine_acceleration" ) = Eigen::Vector3d::Zero( ),
          get_docstring("empirical").c_str());

    m.def("custom", &tss::customAccelerationSettings,
          py::arg( "acceleration_function" ),
          py::arg( "scaling_function" ) = nullptr,
          get_docstring("custom").c_str());

    m.def("direct_tidal_dissipation_acceleration", &tss::directTidalDissipationAcceleration,
          py::arg("k2_love_number"),
          py::arg("time_lag"),
          py::arg("include_direct_radial_component") = true,
          py::arg("use_tide_raised_on_planet") = true,
          get_docstring("direct_tidal_dissipation_acceleration").c_str());

    m.def("momentum_wheel_desaturation_acceleration", &tss::momentumWheelDesaturationAcceleration,
          py::arg("thrust_mid_times"),
          py::arg("delta_v_values"),
          py::arg("total_maneuver_time"),
          py::arg("maneuver_rise_time"),
          get_docstring("momentum_wheel_desaturation_acceleration").c_str());

    m.def("thrust_acceleration", py::overload_cast<const std::shared_ptr<tss::ThrustDirectionSettings>,
                  const std::shared_ptr<tss::ThrustMagnitudeSettings>>(&tss::thrustAcceleration),
          py::arg("thrust_direction_settings"),
          py::arg("thrust_magnitude_settings"),
          get_docstring("thrust_acceleration", 0).c_str());

    m.def("thrust_acceleration", py::overload_cast<
                  const std::shared_ptr<tinterp::DataInterpolationSettings<double, Eigen::Vector3d>>&,
                  const std::function<double(const double)>,
                  const tss::ThrustFrames,
                  const std::string>(&tss::thrustAcceleration),
          py::arg("data_interpolation_settings"),
          py::arg("specific_impulse_function"),
          py::arg("thrust_frame") = tss::ThrustFrames::unspecified_thrust_frame,
          py::arg("central_body") = "",
          get_docstring("thrust_acceleration", 1).c_str());

    m.def("thrust_acceleration", py::overload_cast<
                  const std::shared_ptr<tinterp::DataInterpolationSettings<double, Eigen::Vector3d>>&,
                  const double,
                  const tss::ThrustFrames,
                  const std::string>(&tss::thrustAcceleration),
          py::arg("data_interpolation_settings"),
          py::arg("constant_specific_impulse"),
          py::arg("thrust_frame") = tss::ThrustFrames::unspecified_thrust_frame,
          py::arg("central_body") = "",
          get_docstring("thrust_acceleration", 2).c_str());


    py::class_<tss::AccelerationSettings,
            std::shared_ptr<tss::AccelerationSettings>>(m, "AccelerationSettings",
                    get_docstring("AccelerationSettings").c_str());
//            .def(py::init<const tudat::basic_astrodynamics::AvailableAcceleration>(),
//                 py::arg("acceleration_type"));

    py::class_<tss::SphericalHarmonicAccelerationSettings,
            std::shared_ptr<tss::SphericalHarmonicAccelerationSettings>,
            tss::AccelerationSettings>(m, "SphericalHarmonicAccelerationSettings",
                                       get_docstring("SphericalHarmonicAccelerationSettings").c_str());
//            .def(py::init<const int, const int>(), py::arg("maximum_degree"),
//                 py::arg("maximum_order"));

    py::class_<tss::MutualSphericalHarmonicAccelerationSettings,
            std::shared_ptr<tss::MutualSphericalHarmonicAccelerationSettings>,
            tss::AccelerationSettings>(m, "MutualSphericalHarmonicAccelerationSettings",
                                       get_docstring("MutualSphericalHarmonicAccelerationSettings").c_str());


    py::class_<tss::EmpiricalAccelerationSettings,
            std::shared_ptr<tss::EmpiricalAccelerationSettings>,
            tss::AccelerationSettings>(m, "EmpiricalAccelerationSettings",
                                       get_docstring("EmpiricalAccelerationSettings").c_str());


    py::class_<tss::RelativisticAccelerationCorrectionSettings,
            std::shared_ptr<tss::RelativisticAccelerationCorrectionSettings>,
            tss::AccelerationSettings>(m, "RelativisticAccelerationCorrectionSettings",
                                       get_docstring("RelativisticAccelerationCorrectionSettings").c_str());


    py::class_<tss::CustomAccelerationSettings,
            std::shared_ptr<tss::CustomAccelerationSettings>,
            tss::AccelerationSettings>(m, "CustomAccelerationSettings",
                                       get_docstring("CustomAccelerationSettings").c_str());
    ;

    py::class_<tss::DirectTidalDissipationAccelerationSettings,
            std::shared_ptr<tss::DirectTidalDissipationAccelerationSettings>,
            tss::AccelerationSettings>(m, "DirectTidalDissipationAcceleration",
                                       get_docstring("DirectTidalDissipationAcceleration").c_str());

    py::class_<tss::MomentumWheelDesaturationAccelerationSettings,
            std::shared_ptr<tss::MomentumWheelDesaturationAccelerationSettings>,
            tss::AccelerationSettings>(m, "MomentumWheelDesaturationAcceleration",
                                       get_docstring("MomentumWheelDesaturationAcceleration").c_str());

    py::class_<tss::ThrustAccelerationSettings,
            std::shared_ptr<tss::ThrustAccelerationSettings>,
            tss::AccelerationSettings>(m, "ThrustAccelerationSettings",
                                       get_docstring("ThrustAccelerationSettings").c_str());
//            .def(py::init<//ctor 1
//                         const std::shared_ptr<tss::ThrustDirectionSettings>,
//                         const std::shared_ptr<tss::ThrustMagnitudeSettings>>(),
//                 py::arg("thrust_direction_settings"),
//                 py::arg("thrust_magnitude_settings"))
//            .def(py::init<//ctor 2
//                         const std::shared_ptr<tinterp::DataInterpolationSettings<double, Eigen::Vector3d>> &,
//                         const std::function<double(const double)>,
//                         const tss::ThrustFrames,
//                         const std::string>(),
//                 py::arg("data_interpolation_settings"),
//                 py::arg("specific_impulse_function"),
//                 py::arg("thrust_frame"),
//                 py::arg("central_body") = "")
//            .def(py::init<//ctor 3
//                         const std::shared_ptr<tinterp::DataInterpolationSettings<double, Eigen::Vector3d>> &,
//                         const double,
//                         const tss::ThrustFrames,
//                         const std::string>(),
//                 py::arg("data_interpolation_settings"),
//                 py::arg("constant_specific_impulse"),
//                 py::arg("thrust_frame"),
//                 py::arg("central_body") = "")
            .def_readwrite("direction_settings", &tss::ThrustAccelerationSettings::thrustMagnitudeSettings_ )
            .def_readwrite("magnitude_settings", &tss::ThrustAccelerationSettings::thrustDirectionSettings_ );


    //////////////////////////////////////////////////////////////////////////////
    // createThrustModelGuidance.h / createThrustModelGuidance.cpp
    //////////////////////////////////////////////////////////////////////////////
// For now, it was decided that these functions should not be part of the user interface

//    m.def("get_combined_thrust_direction",
//          &tss::getCombinedThrustDirection,
//          py::arg("thrust_directions"),
//          py::arg("thrust_magnitudes"));
//
//    m.def("get_body_fixed_thrust_direction",
//          &tss::getBodyFixedThrustDirection,
//          py::arg("thrust_magnitude_settings"),
//          py::arg("body_system"),
//          py::arg("body_name"));
//
//    m.def("create_thrust_magnitude_wrapper",
//          &tss::createThrustMagnitudeWrapper,
//          py::arg("thrust_magnitude_settings"),
//          py::arg("body_system"),
//          py::arg("name_of_body_with_guidance"),
//          py::arg("magnitude_update_settings"));
//
//    m.def("update_thrust_settings",
//          &tss::updateThrustSettingsTime,
//          py::arg("thrust_magnitude_wrapper"),
//          py::arg("thrust_direction_guidance"),
//          py::arg("current_time"));
//
//    m.def("reset_thrust_settings_time",
//          &tss::resetThrustSettingsTime,
//          py::arg("thrust_magnitude_wrapper"),
//          py::arg("thrust_direction_guidance"),
//          py::arg("current_time"));

    //////////////////////////////////////////////////////////////////////////////
    // thrustSettings.h / thrustSettings.cpp
    //////////////////////////////////////////////////////////////////////////////
    py::enum_<tss::ThrustDirectionTypes>(m, "ThrustDirectionGuidanceTypes")
            .value("colinear_with_state_segment_thrust_direction", tss::ThrustDirectionTypes::colinear_with_state_segment_thrust_direction)
            .value("thrust_direction_from_existing_body_orientation", tss::ThrustDirectionTypes::thrust_direction_from_existing_body_orientation)
            .value("custom_thrust_direction", tss::ThrustDirectionTypes::custom_thrust_direction)
            .value("custom_thrust_orientation", tss::ThrustDirectionTypes::custom_thrust_orientation)
            .value("mee_costate_based_thrust_direction", tss::ThrustDirectionTypes::mee_costate_based_thrust_direction);

    m.def("get_propulsion_input_variables",
          &tss::getPropulsionInputVariables,
          py::arg("body_with_guidance") = std::shared_ptr<tss::Body>(),
          py::arg("independent_variables") = std::vector<tudat::propulsion::ThrustIndependentVariables>(),
          py::arg("guidance_input_functions") = std::vector<std::function<double()>>(),
          get_docstring("get_propulsion_input_variables").c_str());

    py::class_<
            tss::ThrustDirectionSettings,
            std::shared_ptr<tss::ThrustDirectionSettings>>(m, "ThrustDirectionSettings",
                    get_docstring("ThrustDirectionSettings").c_str())
//            .def(py::init<
//                         const tss::ThrustDirectionGuidanceTypes,
//                         const std::string>(),
//                 py::arg("thrust_direction_type"),
//                 py::arg("relative_body"))
            .def_readonly("thrust_direction_type", &tss::ThrustDirectionSettings::thrustDirectionType_)
            .def_readonly("relative_body", &tss::ThrustDirectionSettings::relativeBody_);

    py::class_<
            tss::ThrustMagnitudeSettings,
            std::shared_ptr<tss::ThrustMagnitudeSettings>>(m, "ThrustMagnitudeSettings",
                                                           get_docstring("ThrustMagnitudeSettings").c_str())
//            .def(py::init<
//                         const tss::ThrustMagnitudeTypes,
//                         const std::string &>(),
//                 py::arg("thrust_magnitude_guidance_type"),
//                 py::arg("thrust_origin_id"))
            .def_readonly("thrust_magnitude_type", &tss::ThrustMagnitudeSettings::thrustMagnitudeType_)
            .def_readonly("thrust_origin_id", &tss::ThrustMagnitudeSettings::thrustOriginId_);

    py::class_<
            tss::ThrustDirectionFromStateGuidanceSettings,
            std::shared_ptr<tss::ThrustDirectionFromStateGuidanceSettings>,
            tss::ThrustDirectionSettings>(m, "ThrustDirectionFromStateGuidanceSettings",
                                          get_docstring("ThrustDirectionFromStateGuidanceSettings").c_str())
//            .def(py::init<const std::string &,
//                         const bool,
//                         const bool>(),
//                 py::arg("central_body"),
//                 py::arg("is_colinear_with_velocity"),
//                 py::arg("direction_is_opposite_to_vector"))
            .def_readonly("is_colinear_with_velocity", &tss::ThrustDirectionFromStateGuidanceSettings::isColinearWithVelocity_)
            .def_readonly("direction_is_opposite_to_vector", &tss::ThrustDirectionFromStateGuidanceSettings::directionIsOppositeToVector_);

    py::class_<
            tss::CustomThrustDirectionSettings,
            std::shared_ptr<tss::CustomThrustDirectionSettings>,
            tss::ThrustDirectionSettings>(m, "CustomThrustDirectionSettings",
                                          get_docstring("CustomThrustDirectionSettings").c_str())
//            .def(py::init<const std::function<Eigen::Vector3d(const double)>>(),
//                 py::arg("thrust_direction_function") );
            .def_readonly("thrust_direction_function", &tss::CustomThrustDirectionSettings::thrustDirectionFunction_);

    py::class_<
            tss::CustomThrustOrientationSettings,
            std::shared_ptr<tss::CustomThrustOrientationSettings>,
            tss::ThrustDirectionSettings>(m, "CustomThrustOrientationSettings",
                                          get_docstring("CustomThrustOrientationSettings").c_str())
//            .def(py::init<const std::function<Eigen::Matrix3d(const double)>>(),
//                 py::arg("thrust_orientation_function"))
            .def_readonly("thrust_orientation_function", &tss::CustomThrustOrientationSettings::thrustOrientationFunction_);

    py::class_<
            tss::MeeCostateBasedThrustDirectionSettings,
            std::shared_ptr<tss::MeeCostateBasedThrustDirectionSettings>,
            tss::ThrustDirectionSettings>(m, "MeeCostateBasedThrustDirectionSettings",
                                          get_docstring("MeeCostateBasedThrustDirectionSettings").c_str())
//            .def(py::init<const std::string &,//ctor 1
//                         const std::string &,
//                         const std::function<Eigen::VectorXd(const double)>>(),
//                 py::arg("vehicle_name"),
//                 py::arg("central_body_name"),
//                 py::arg("costate_function"))
//            .def(py::init<const std::string &,//ctor 2
//                         const std::string &,
//                         std::shared_ptr<tinterp::OneDimensionalInterpolator<double, Eigen::VectorXd>>>(),
//                 py::arg("vehicle_name"),
//                 py::arg("central_body_name"),
//                 py::arg("costate_interpolator"))
//            .def(py::init<const std::string &,//ctor 3
//                         const std::string &,
//                         const Eigen::VectorXd>(),
//                 py::arg("vehicle_name"),
//                 py::arg("central_body_name"),
//                 py::arg("constant_costates"))
            .def_readonly("vehicle_name", &tss::MeeCostateBasedThrustDirectionSettings::vehicleName_)
            .def_readonly("costate_function", &tss::MeeCostateBasedThrustDirectionSettings::costateFunction_);

    py::enum_<tss::ThrustMagnitudeTypes>(m, "ThrustMagnitudeTypes")
            .value("constant_thrust_magnitude", tss::ThrustMagnitudeTypes::constant_thrust_magnitude)
            .value("from_engine_properties_thrust_magnitude", tss::ThrustMagnitudeTypes::from_engine_properties_thrust_magnitude)
            .value("thrust_magnitude_from_time_function", tss::ThrustMagnitudeTypes::thrust_magnitude_from_time_function)
            .value("thrust_magnitude_from_dependent_variables", tss::ThrustMagnitudeTypes::thrust_magnitude_from_dependent_variables)
            .value("bang_bang_thrust_magnitude_from_mee_costates", tss::ThrustMagnitudeTypes::bang_bang_thrust_magnitude_from_mee_costates);

    py::class_<
            tss::ConstantThrustMagnitudeSettings,
            std::shared_ptr<tss::ConstantThrustMagnitudeSettings>,
            tss::ThrustMagnitudeSettings>(m, "ConstantThrustMagnitudeSettings",
                                          get_docstring("ConstantThrustMagnitudeSettings").c_str())
//            .def(py::init<
//                         const double,
//                         const double,
//                         const Eigen::Vector3d>(),
//                 py::arg("thrust_magnitude"),
//                 py::arg("specific_impulse"),
//                 py::arg("body_fixed_thrust_direction") = Eigen::Vector3d::UnitX())
            .def_readonly("thrust_magnitude", &tss::ConstantThrustMagnitudeSettings::thrustMagnitude_)
            .def_readonly("specific_impulse", &tss::ConstantThrustMagnitudeSettings::specificImpulse_)
            .def_readonly("body_fixed_thrust_direction", &tss::ConstantThrustMagnitudeSettings::bodyFixedThrustDirection_);

//    py::class_<
//            tss::FromBodyThrustMagnitudeSettings,
//            std::shared_ptr<tss::FromBodyThrustMagnitudeSettings>,
//            tss::ThrustMagnitudeSettings>(m, "FromBodyThrustMagnitudeSettings")
//            .def(py::init<
//                         const double,
//                         const std::string &>(),
//                 py::arg("use_all_engines"),
//                 py::arg("thrust_origin"))
//            .def_readonly("use_all_engines", &tss::FromBodyThrustMagnitudeSettings::useAllEngines_);

    py::class_<
            tss::FromFunctionThrustMagnitudeSettings,
            std::shared_ptr<tss::FromFunctionThrustMagnitudeSettings>,
            tss::ThrustMagnitudeSettings>(m, "FromFunctionThrustMagnitudeSettings",
                                          get_docstring("FromFunctionThrustMagnitudeSettings").c_str())
//            .def(py::init<
//                         const std::function< double( const double ) >,
//                         const std::function< double( const double ) >,
//                         const std::function< bool( const double ) >,
//                         const std::function< Eigen::Vector3d( ) >,
//                         const std::function< void( const double ) > >(),
//                 py::arg("thrust_magnitude_function"),
//                 py::arg("specific_impulse_function"),
//                 py::arg("is_engine_on_function" ) =
//                         std::function< bool( const double ) >( [ ]( const double ){ return true; } ),
//                 py::arg("body_fixed_thrust_direction" ) =
//                         std::function< Eigen::Vector3d( ) >( [ ]( ){ return  Eigen::Vector3d::UnitX( ); } ),
//                 py::arg("custom_thrust_reset_function" ) = std::function< void( const double ) >( ) );

    // Thrust orientation factory functions

    m.def("thrust_direction_from_state_guidance", &tss::thrustDirectionFromStateGuidanceSettings,
          py::arg( "central_body"),
          py::arg("is_colinear_with_velocity"),
          py::arg("direction_is_opposite_to_vector"),
          get_docstring("thrust_direction_from_state_guidance").c_str());

    m.def("thrust_from_existing_body_orientation", &tss::thrustFromExistingBodyOrientation,
          get_docstring("thrust_from_existing_body_orientation").c_str());

    m.def("custom_thrust_orientation",
          py::overload_cast< std::function< Eigen::Matrix3d( const double ) > >(
                  &tss::customThrustOrientationSettings ),
          py::arg( "thrust_orientation_function" ),
          get_docstring("custom_thrust_orientation").c_str());

    m.def("custom_thrust_direction", &tss::customThrustDirectionSettings,
          py::arg( "thrust_direction_function" ),
          get_docstring("custom_thrust_direction").c_str());

    m.def("mee_costate_based_thrust_direction",
          py::overload_cast<const std::string&, const std::string&,
                  const std::shared_ptr< tinterp::OneDimensionalInterpolator< double, Eigen::VectorXd> > >(
                  &tss::meeCostateBasedThrustDirectionSettings),
          py::arg("vehicle_name"),
          py::arg("central_body_name"),
          py::arg("costate_interpolator"),
          get_docstring("mee_costate_based_thrust_direction").c_str());

    m.def("mee_costate_based_thrust_direction",
          py::overload_cast<const std::string&, const std::string&,
                  const Eigen::VectorXd >(
                  &tss::meeCostateBasedThrustDirectionSettings),
          py::arg("vehicle_name"),
          py::arg("central_body_name"),
          py::arg("constant_costates"),
          get_docstring("mee_costate_based_thrust_direction").c_str());

    // Thrust orientation factory functions

    m.def("constant_thrust_magnitude", &tss::constantThrustMagnitudeSettings,
          py::arg("thrust_magnitude"),
          py::arg("specific_impulse"),
          py::arg("body_fixed_thrust_direction" ) =
                  Eigen::Vector3d::UnitX(),
          get_docstring("constant_thrust_magnitude").c_str());

    m.def("custom_thrust_magnitude", &tss::fromFunctionThrustMagnitudeSettings,
          py::arg("thrust_magnitude_function"),
          py::arg("specific_impulse_function"),
          py::arg("is_engine_on_function" ) =
                  std::function< bool( const double ) >( [ ]( const double ){ return true; } ),
          py::arg("body_fixed_thrust_direction" ) =
                  std::function< Eigen::Vector3d( ) >( [ ]( ){ return  Eigen::Vector3d::UnitX( ); } ),
          py::arg("custom_thrust_reset_function" ) = std::function< void( const double ) >( ),
          get_docstring("custom_thrust_magnitude").c_str());

    // TODO: EngineModel still to be implemented
//    m.def("from_body_thrust_magnitude", &tss::fromBodyThrustMagnitudeSettings,
//          py::arg("use_all_engines") = false,
//          py::arg("thrust_origin") = "");

}

}// namespace propagation_setup
}// namespace simulation
}// namespace tudatpy
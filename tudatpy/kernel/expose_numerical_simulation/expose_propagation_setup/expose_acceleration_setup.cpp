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
//#include "kernel/expose_numerical_simulation/deprecation_support.h"

#include <tudat/basics/deprecationWarnings.h>

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

namespace tudat {
namespace simulation_setup {


inline std::shared_ptr< AccelerationSettings > customAccelerationSettingsDeprecated(
        const std::function< Eigen::Vector3d( const double ) > accelerationFunction )
{
    static bool isWarningPrinted = false;
    if( isWarningPrinted == false )
    {
        tudat::utilities::printDeprecationWarning( "tudatpy.numerical_simulation.propagation_setup.acceleration.custom",
                             "tudatpy.numerical_simulation.propagation_setup.acceleration.custom_acceleration");
        isWarningPrinted = true;
    }

    return customAccelerationSettings( accelerationFunction );
}


inline std::shared_ptr< AccelerationSettings > thrustAccelerationRemoved1(
        const std::shared_ptr<tss::ThrustDirectionSettings> thrustDirectionSettings,
        const std::shared_ptr<tss::ThrustMagnitudeSettings> thrustMagnitudeSettings )
{
    tudat::utilities::printDeprecationError( "tudatpy.numerical_simulation.propagation_setup.acceleration.thrust_from_direction_and_magnitude",
                                             "https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/environment_setup/thrust_refactor/thrust_refactor.html#thrust-acceleration" );
    return nullptr;
}

inline std::shared_ptr< AccelerationSettings > thrustAccelerationRemoved2(
        const std::function< Eigen::Vector3d( const double ) > thrustForceFunction,
        const std::function< double( const double ) > specificImpulseFunction,
        const ThrustFrames thrustFrame = unspecified_thrust_frame,
        const std::string centralBody = "" )
{
    tudat::utilities::printDeprecationError( "tudatpy.numerical_simulation.propagation_setup.acceleration.thrust_and_isp_from_custom_function",
                                             "https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/environment_setup/thrust_refactor/thrust_refactor.html#thrust-acceleration" );
    return nullptr;

}

inline std::shared_ptr< AccelerationSettings > thrustAccelerationRemoved3(
        const std::function< Eigen::Vector3d( const double ) > thrustForceFunction,
        const double constantSpecificImpulse,
        const ThrustFrames thrustFrame = unspecified_thrust_frame,
        const std::string centralBody = "" )
{
    tudat::utilities::printDeprecationError( "tudatpy.numerical_simulation.propagation_setup.acceleration.thrust_from_custom_function",
                                             "https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/environment_setup/thrust_refactor/thrust_refactor.html#thrust-acceleration" );
    return nullptr;
}


//! @get_docstring(customAccelerationSettings)
inline std::shared_ptr< AccelerationSettings > customAccelerationSettings(
        const std::function< Eigen::Vector3d( const double ) > accelerationFunction )
{
        return std::make_shared< CustomAccelerationSettings >(
                    accelerationFunction );

}

}

}

namespace tudatpy {
namespace numerical_simulation {
namespace propagation_setup {
namespace acceleration {

void expose_acceleration_setup(py::module &m) {

    /*
     * This contains the addition of IntegratorSettings and AvailableIntegrators
     * and AvailableAccelerations which should be relocated in the tudat source.
     */


    py::enum_<tba::AvailableAcceleration>(m, "AvailableAcceleration", get_docstring("AvailableAcceleration").c_str())
            .value("undefined_acceleration_type", tba::AvailableAcceleration::undefined_acceleration, get_docstring("AvailableAcceleration.undefined_acceleration_type").c_str())
            .value("point_mass_gravity_type", tba::AvailableAcceleration::point_mass_gravity, get_docstring("AvailableAcceleration.point_mass_gravity_type").c_str())
            .value("aerodynamic_type", tba::AvailableAcceleration::aerodynamic, get_docstring("AvailableAcceleration.aerodynamic_type").c_str())
            .value("cannonball_radiation_pressure_type", tba::AvailableAcceleration::cannon_ball_radiation_pressure, get_docstring("AvailableAcceleration.cannon_ball_radiation_pressure_type").c_str())
            .value("spherical_harmonic_gravity_type", tba::AvailableAcceleration::spherical_harmonic_gravity, get_docstring("AvailableAcceleration.spherical_harmonic_gravity_type").c_str())
            .value("mutual_spherical_harmonic_gravity_type", tba::AvailableAcceleration::mutual_spherical_harmonic_gravity, get_docstring("AvailableAcceleration.mutual_spherical_harmonic_gravity_type").c_str())
            .value("polyhedron_gravity_type", tba::AvailableAcceleration::polyhedron_gravity, get_docstring("AvailableAcceleration.polyhedron_gravity_type").c_str())
            .value("ring_gravity_type", tba::AvailableAcceleration::ring_gravity, get_docstring("AvailableAcceleration.ring_gravity_type").c_str())
            .value("thrust_acceleration_type", tba::AvailableAcceleration::thrust_acceleration, get_docstring("AvailableAcceleration.thrust_acceleration_type").c_str())
            .value("relativistic_correction_acceleration_type", tba::AvailableAcceleration::relativistic_correction_acceleration, get_docstring("AvailableAcceleration.relativistic_correction_acceleration_type").c_str())
            .value("empirical_acceleration_type", tba::AvailableAcceleration::empirical_acceleration, get_docstring("AvailableAcceleration.empirical_acceleration_type").c_str())
            .value("direct_tidal_dissipation_in_central_body_acceleration_type", tba::AvailableAcceleration::direct_tidal_dissipation_in_central_body_acceleration, get_docstring("AvailableAcceleration.direct_tidal_dissipation_in_central_body_acceleration_type").c_str())
            .value("direct_tidal_dissipation_in_orbiting_body_acceleration_type", tba::AvailableAcceleration::direct_tidal_dissipation_in_orbiting_body_acceleration, get_docstring("AvailableAcceleration.direct_tidal_dissipation_in_orbiting_body_acceleration_type").c_str())
            .value("quasi_impulsive_shots_acceleration_type", tba::AvailableAcceleration::momentum_wheel_desaturation_acceleration, get_docstring("AvailableAcceleration.quasi_impulsive_shots_acceleration_type").c_str())
            .value("custom_acceleration_type", tba::AvailableAcceleration::custom_acceleration, get_docstring("AvailableAcceleration.custom_acceleration_type").c_str())
            .value("radiation_pressure_type", tba::AvailableAcceleration::radiation_pressure, get_docstring("AvailableAcceleration.radiation_pressure_type").c_str())
            .export_values();

    //////////////////////////////////////////////////////////////////////////////
    // accelerationSettings.h
    //////////////////////////////////////////////////////////////////////////////

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


    py::class_<tss::DirectTidalDissipationAccelerationSettings,
            std::shared_ptr<tss::DirectTidalDissipationAccelerationSettings>,
            tss::AccelerationSettings>(m, "DirectTidalDissipationAccelerationSettings",
                                       get_docstring("DirectTidalDissipationAccelerationSettings").c_str());

    py::class_<tss::MomentumWheelDesaturationAccelerationSettings,
            std::shared_ptr<tss::MomentumWheelDesaturationAccelerationSettings>,
            tss::AccelerationSettings>(m, "MomentumWheelDesaturationAccelerationSettings",
                                       get_docstring("MomentumWheelDesaturationAccelerationSettings").c_str());

    py::class_<tss::ThrustAccelerationSettings,
            std::shared_ptr<tss::ThrustAccelerationSettings>,
            tss::AccelerationSettings>(m, "ThrustAccelerationSettings",
                                       get_docstring("ThrustAccelerationSettings").c_str())
            .def_property_readonly("direction_settings",
                                   &tss::ThrustAccelerationSettings::printDeprecationError<
                                   std::shared_ptr< tss::ThrustDirectionSettings > > )
            .def_property_readonly("magnitude_settings",
                                   &tss::ThrustAccelerationSettings::printDeprecationError<
                                   std::shared_ptr< tss::ThrustMagnitudeSettings > > );


    // Unified interface functions for acceleration settings
    //  m.def("acceleration", &tss::acceleration, py::arg("acceleration_type"));
    m.def("point_mass_gravity", &tss::pointMassGravityAcceleration,
          get_docstring("point_mass_gravity").c_str());

    m.def("aerodynamic", &tss::aerodynamicAcceleration,
          get_docstring("aerodynamic").c_str());

    m.def("cannonball_radiation_pressure", &tss::cannonBallRadiationPressureAcceleration,
          get_docstring("cannonball_radiation_pressure").c_str());

    m.def("radiation_pressure", &tss::radiationPressureAcceleration,
          get_docstring("radiation_pressure").c_str());

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

     m.def("polyhedron_gravity", &tss::polyhedronAcceleration,
          get_docstring("polyhedron_gravity").c_str());

    m.def("ring_gravity", &tss::ringAcceleration,
          get_docstring("ring_gravity").c_str());

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

    m.def("yarkovsky", &tss::yarkovskyAcceleration,
          py::arg( "yarkovsky_parameter" ),
          get_docstring("yarkovsky").c_str());

    m.def("custom",
          &tss::customAccelerationSettingsDeprecated,
          py::arg( "acceleration_function" ) );


    m.def("custom_acceleration",
          py::overload_cast< std::function< Eigen::Vector3d( const double ) > >(
              &tss::customAccelerationSettings ),
          py::arg( "acceleration_function" ),
          get_docstring("custom_acceleration").c_str());

    m.def("direct_tidal_dissipation_acceleration", &tss::directTidalDissipationAcceleration,
          py::arg("k2_love_number"),
          py::arg("time_lag"),
          py::arg("include_direct_radial_component") = true,
          py::arg("use_tide_raised_on_planet") = true,
          py::arg("explicit_libraional_tide_on_satellite" ) = false,
          get_docstring("direct_tidal_dissipation_acceleration").c_str());

    m.def("direct_tidal_dissipation_acceleration", &tss::directTidalDissipationAccelerationFromInvQ,
          py::arg("k2_love_number"),
          py::arg("inverse_tidal_quality_factor"),
          py::arg("tidal_period"),
          py::arg("include_direct_radial_component") = true,
          py::arg("use_tide_raised_on_planet") = true,
          py::arg("explicit_libraional_tide_on_satellite" ) = false,
          get_docstring("direct_tidal_dissipation_acceleration").c_str());

    m.def("quasi_impulsive_shots_acceleration", &tss::momentumWheelDesaturationAcceleration,
          py::arg("thrust_mid_times"),
          py::arg("delta_v_values"),
          py::arg("total_maneuver_time"),
          py::arg("maneuver_rise_time"),
          get_docstring("quasi_impulsive_shots_acceleration").c_str());

    m.def("thrust_from_engines", &tss::thrustAcceleration,
          py::arg("engine_names"),
          get_docstring("thrust_from_engines").c_str());

    m.def("thrust_from_engine", &tss::thrustAccelerationFromSingleEngine,
          py::arg("engine_name"),
          get_docstring("thrust_from_engine").c_str());


    m.def("thrust_from_all_engines", &tss::thrustAccelerationFromAllEngines,
          get_docstring("thrust_from_all_engines").c_str());



    m.def("thrust_from_direction_and_magnitude", &tss::thrustAccelerationRemoved1,
          py::arg("thrust_direction_settings"),
          py::arg("thrust_magnitude_settings") );

    m.def("thrust_from_custom_function", &tss::thrustAccelerationRemoved2,
          py::arg("thrust_force_function"),
          py::arg("specific_impulse_function"),
          py::arg("thrust_frame") = tss::ThrustFrames::inertial_thrust_frame,
          py::arg("central_body") = "" );

    m.def("thrust_and_isp_from_custom_function", &tss::thrustAccelerationRemoved3,
          py::arg("thrust_force_function"),
          py::arg("constant_specific_impulse"),
          py::arg("thrust_frame") = tss::ThrustFrames::inertial_thrust_frame,
          py::arg("central_body") = "",
          get_docstring("thrust_and_isp_from_custom_function").c_str());


}

}// namespace acceleration
}// namespace propagation_setup
}// namespace numerical_simulation
}// namespace tudatpy

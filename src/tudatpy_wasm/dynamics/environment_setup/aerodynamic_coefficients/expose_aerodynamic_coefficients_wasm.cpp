/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include "../../../wasm_module.h"
#include "../../../eigen_wasm.h"
#include "../../../shared_ptr_wasm.h"

#include <tudat/simulation/environment_setup/createAerodynamicCoefficientInterface.h>
#include <tudat/astro/aerodynamics/aerodynamicCoefficientInterface.h>
#include <tudat/astro/reference_frames/referenceFrameTransformations.h>

namespace tss = tudat::simulation_setup;
namespace ta = tudat::aerodynamics;
namespace trf = tudat::reference_frames;

WASM_MODULE_PATH("dynamics_environment_setup_aerodynamic_coefficients")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_environment_setup_aerodynamic_coefficients) {
    using namespace emscripten;

    // AerodynamicCoefficientsIndependentVariables enum
    enum_<ta::AerodynamicCoefficientsIndependentVariables>(
        "dynamics_environment_setup_aerodynamic_coefficients_AerodynamicCoefficientsIndependentVariables")
        .value("mach_number_dependent", ta::mach_number_dependent)
        .value("angle_of_attack_dependent", ta::angle_of_attack_dependent)
        .value("sideslip_angle_dependent", ta::angle_of_sideslip_dependent)
        .value("altitude_dependent", ta::altitude_dependent)
        .value("time_dependent", ta::time_dependent)
        .value("control_surface_deflection_dependent", ta::control_surface_deflection_dependent)
        .value("undefined_independent_variable", ta::undefined_independent_variable);

    // AerodynamicCoefficientFrames enum
    enum_<ta::AerodynamicCoefficientFrames>(
        "dynamics_environment_setup_aerodynamic_coefficients_AerodynamicCoefficientFrames")
        .value("positive_body_fixed_frame_coefficients", ta::body_fixed_frame_coefficients)
        .value("negative_body_fixed_frame_coefficients", ta::negative_body_fixed_frame_coefficients)
        .value("positive_aerodynamic_frame_coefficients", ta::positive_aerodynamic_frame_coefficients)
        .value("negative_aerodynamic_frame_coefficients", ta::negative_aerodynamic_frame_coefficients);

    // AerodynamicsReferenceFrames enum
    enum_<trf::AerodynamicsReferenceFrames>(
        "dynamics_environment_setup_aerodynamic_coefficients_AerodynamicsReferenceFrames")
        .value("inertial_frame", trf::inertial_frame)
        .value("corotating_frame", trf::corotating_frame)
        .value("vertical_frame", trf::vertical_frame)
        .value("trajectory_frame", trf::trajectory_frame)
        .value("aerodynamic_frame", trf::aerodynamic_frame)
        .value("body_frame", trf::body_frame);

    // AerodynamicsReferenceFrameAngles enum
    enum_<trf::AerodynamicsReferenceFrameAngles>(
        "dynamics_environment_setup_aerodynamic_coefficients_AerodynamicsReferenceFrameAngles")
        .value("latitude_angle", trf::latitude_angle)
        .value("longitude_angle", trf::longitude_angle)
        .value("heading_angle", trf::heading_angle)
        .value("flight_path_angle", trf::flight_path_angle)
        .value("angle_of_attack", trf::angle_of_attack)
        .value("angle_of_sideslip", trf::angle_of_sideslip)
        .value("bank_angle", trf::bank_angle);

    // AerodynamicCoefficientSettings base class
    class_<tss::AerodynamicCoefficientSettings>(
        "dynamics_environment_setup_aerodynamic_coefficients_AerodynamicCoefficientSettings")
        .smart_ptr<std::shared_ptr<tss::AerodynamicCoefficientSettings>>(
            "shared_ptr_AerodynamicCoefficientSettings");

    // ConstantAerodynamicCoefficientSettings derived class
    class_<tss::ConstantAerodynamicCoefficientSettings, base<tss::AerodynamicCoefficientSettings>>(
        "dynamics_environment_setup_aerodynamic_coefficients_ConstantAerodynamicCoefficientSettings")
        .smart_ptr<std::shared_ptr<tss::ConstantAerodynamicCoefficientSettings>>(
            "shared_ptr_ConstantAerodynamicCoefficientSettings");

    // Factory functions
    function("dynamics_environment_setup_aerodynamic_coefficients_constant",
        select_overload<std::shared_ptr<tss::AerodynamicCoefficientSettings>(
            const double,
            const Eigen::Vector3d&,
            const ta::AerodynamicCoefficientFrames)>(&tss::constantAerodynamicCoefficientSettings));

    function("dynamics_environment_setup_aerodynamic_coefficients_constant_force_and_moment",
        &tss::constantAerodynamicForceAndMomentCoefficientSettings);
}

#endif

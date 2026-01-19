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
#include "../../../stl_wasm.h"
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

    // ========================================================================
    // Enumerations
    // ========================================================================

    // AerodynamicCoefficientsIndependentVariables enum
    enum_<ta::AerodynamicCoefficientsIndependentVariables>(
        "dynamics_environment_setup_aerodynamic_coefficients_AerodynamicCoefficientsIndependentVariables")
        .value("mach_number_dependent", ta::mach_number_dependent)
        .value("angle_of_attack_dependent", ta::angle_of_attack_dependent)
        .value("sideslip_angle_dependent", ta::angle_of_sideslip_dependent)
        .value("altitude_dependent", ta::altitude_dependent)
        .value("time_dependent", ta::time_dependent)
        .value("temperature_dependent", ta::temperature_dependent)
        .value("velocity_dependent", ta::velocity_dependent)
        .value("he_number_density_dependent", ta::he_number_density_dependent)
        .value("o_number_density_dependent", ta::o_number_density_dependent)
        .value("n2_number_density_dependent", ta::n2_number_density_dependent)
        .value("o2_number_density_dependent", ta::o2_number_density_dependent)
        .value("ar_number_density_dependent", ta::ar_number_density_dependent)
        .value("h_number_density_dependent", ta::h_number_density_dependent)
        .value("n_number_density_dependent", ta::n_number_density_dependent)
        .value("anomalous_o_number_density_dependent", ta::anomalous_o_number_density_dependent)
        .value("control_surface_deflection_dependent", ta::control_surface_deflection_dependent)
        .value("undefined_independent_variable", ta::undefined_independent_variable);

    // AerodynamicCoefficientFrames enum
    enum_<ta::AerodynamicCoefficientFrames>(
        "dynamics_environment_setup_aerodynamic_coefficients_AerodynamicCoefficientFrames")
        .value("positive_body_fixed_frame_coefficients", ta::body_fixed_frame_coefficients)
        .value("negative_body_fixed_frame_coefficients", ta::negative_body_fixed_frame_coefficients)
        .value("positive_aerodynamic_frame_coefficients", ta::positive_aerodynamic_frame_coefficients)
        .value("negative_aerodynamic_frame_coefficients", ta::negative_aerodynamic_frame_coefficients);

    // AtmosphericCompositionSpecies enum
    enum_<ta::AtmosphericCompositionSpecies>(
        "dynamics_environment_setup_aerodynamic_coefficients_AtmosphericCompositionSpecies")
        .value("o_species", ta::o_species)
        .value("o2_species", ta::o2_species)
        .value("n2_species", ta::n2_species)
        .value("he_species", ta::he_species)
        .value("h_species", ta::h_species)
        .value("ar_species", ta::ar_species)
        .value("n_species", ta::n_species)
        .value("anomalous_o_species", ta::anomalous_o_species);

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

    // GasSurfaceInteractionModelType enum
    enum_<ta::GasSurfaceInteractionModelType>(
        "dynamics_environment_setup_aerodynamic_coefficients_GasSurfaceInteractionModelType")
        .value("newton", ta::newton)
        .value("storch", ta::storch)
        .value("sentman", ta::sentman)
        .value("cook", ta::cook);

    // ========================================================================
    // Settings classes
    // ========================================================================

    // AerodynamicCoefficientSettings base class
    class_<tss::AerodynamicCoefficientSettings>(
        "dynamics_environment_setup_aerodynamic_coefficients_AerodynamicCoefficientSettings")
        .smart_ptr<std::shared_ptr<tss::AerodynamicCoefficientSettings>>(
            "shared_ptr_AerodynamicCoefficientSettings")
        .function("getAddForceContributionToMoments",
            &tss::AerodynamicCoefficientSettings::getAddForceContributionToMoments)
        .function("setAddForceContributionToMoments",
            &tss::AerodynamicCoefficientSettings::setAddForceContributionToMoments)
        .function("getMomentReferencePoint",
            &tss::AerodynamicCoefficientSettings::getMomentReferencePoint)
        .function("setMomentReferencePoint",
            &tss::AerodynamicCoefficientSettings::setMomentReferencePoint)
        .function("addControlSurfaceSettings",
            &tss::AerodynamicCoefficientSettings::addControlSurfaceSettings);

    // ConstantAerodynamicCoefficientSettings
    class_<tss::ConstantAerodynamicCoefficientSettings, base<tss::AerodynamicCoefficientSettings>>(
        "dynamics_environment_setup_aerodynamic_coefficients_ConstantAerodynamicCoefficientSettings")
        .smart_ptr<std::shared_ptr<tss::ConstantAerodynamicCoefficientSettings>>(
            "shared_ptr_ConstantAerodynamicCoefficientSettings");

    // CustomAerodynamicCoefficientSettings
    class_<tss::CustomAerodynamicCoefficientSettings, base<tss::AerodynamicCoefficientSettings>>(
        "dynamics_environment_setup_aerodynamic_coefficients_CustomAerodynamicCoefficientSettings")
        .smart_ptr<std::shared_ptr<tss::CustomAerodynamicCoefficientSettings>>(
            "shared_ptr_CustomAerodynamicCoefficientSettings");

    // ScaledAerodynamicCoefficientInterfaceSettings
    class_<tss::ScaledAerodynamicCoefficientInterfaceSettings, base<tss::AerodynamicCoefficientSettings>>(
        "dynamics_environment_setup_aerodynamic_coefficients_ScaledAerodynamicCoefficientInterfaceSettings")
        .smart_ptr<std::shared_ptr<tss::ScaledAerodynamicCoefficientInterfaceSettings>>(
            "shared_ptr_ScaledAerodynamicCoefficientInterfaceSettings");

    // ControlSurfaceIncrementAerodynamicCoefficientSettings
    class_<tss::ControlSurfaceIncrementAerodynamicCoefficientSettings>(
        "dynamics_environment_setup_aerodynamic_coefficients_ControlSurfaceIncrementAerodynamicCoefficientSettings")
        .smart_ptr<std::shared_ptr<tss::ControlSurfaceIncrementAerodynamicCoefficientSettings>>(
            "shared_ptr_ControlSurfaceIncrementAerodynamicCoefficientSettings");

    // CustomControlSurfaceIncrementAerodynamicCoefficientSettings
    class_<tss::CustomControlSurfaceIncrementAerodynamicCoefficientSettings,
           base<tss::ControlSurfaceIncrementAerodynamicCoefficientSettings>>(
        "dynamics_environment_setup_aerodynamic_coefficients_CustomControlSurfaceIncrementAerodynamicCoefficientSettings")
        .smart_ptr<std::shared_ptr<tss::CustomControlSurfaceIncrementAerodynamicCoefficientSettings>>(
            "shared_ptr_CustomControlSurfaceIncrementAerodynamicCoefficientSettings");

    // ========================================================================
    // Factory functions
    // ========================================================================

    // Constant coefficients (force only)
    function("dynamics_environment_setup_aerodynamic_coefficients_constant",
        select_overload<std::shared_ptr<tss::AerodynamicCoefficientSettings>(
            const double,
            const Eigen::Vector3d&,
            const ta::AerodynamicCoefficientFrames)>(&tss::constantAerodynamicCoefficientSettings));

    // Constant coefficients (force and moment)
    function("dynamics_environment_setup_aerodynamic_coefficients_constant_force_and_moment",
        &tss::constantAerodynamicForceAndMomentCoefficientSettings);

    // Custom force coefficients (force only)
    function("dynamics_environment_setup_aerodynamic_coefficients_custom_aerodynamic_force_coefficients",
        select_overload<std::shared_ptr<tss::AerodynamicCoefficientSettings>(
            const std::function<Eigen::Vector3d(const std::vector<double>&)>,
            const double,
            const std::vector<ta::AerodynamicCoefficientsIndependentVariables>,
            const ta::AerodynamicCoefficientFrames)>(&tss::customAerodynamicCoefficientSettings));

    // Custom force and moment coefficients
    function("dynamics_environment_setup_aerodynamic_coefficients_custom_aerodynamic_force_and_moment_coefficients",
        select_overload<std::shared_ptr<tss::AerodynamicCoefficientSettings>(
            const std::function<Eigen::Vector3d(const std::vector<double>&)>,
            const std::function<Eigen::Vector3d(const std::vector<double>&)>,
            const double,
            const double,
            const std::vector<ta::AerodynamicCoefficientsIndependentVariables>,
            const ta::AerodynamicCoefficientFrames,
            const ta::AerodynamicCoefficientFrames,
            const Eigen::Vector3d&)>(&tss::customAerodynamicCoefficientSettings));

    // Scaled coefficients by constant
    function("dynamics_environment_setup_aerodynamic_coefficients_scaled_by_constant",
        select_overload<std::shared_ptr<tss::AerodynamicCoefficientSettings>(
            const std::shared_ptr<tss::AerodynamicCoefficientSettings>,
            const double,
            const double,
            const bool)>(&tss::scaledAerodynamicCoefficientSettings));

    // Scaled coefficients by vector
    function("dynamics_environment_setup_aerodynamic_coefficients_scaled_by_vector",
        select_overload<std::shared_ptr<tss::AerodynamicCoefficientSettings>(
            const std::shared_ptr<tss::AerodynamicCoefficientSettings>,
            const Eigen::Vector3d,
            const Eigen::Vector3d,
            const bool)>(&tss::scaledAerodynamicCoefficientSettings));

    // Panelled aerodynamic coefficients
    function("dynamics_environment_setup_aerodynamic_coefficients_panelled",
        select_overload<std::shared_ptr<tss::AerodynamicCoefficientSettings>(
            const ta::GasSurfaceInteractionModelType,
            const double,
            const int,
            const bool)>(&tss::panelledAerodynamicCoefficientSettings));

    // Custom control surface increment coefficients
    function("dynamics_environment_setup_aerodynamic_coefficients_custom_control_surface",
        &tss::customControlSurfaceIncrementAerodynamicCoefficientSettings);
}

#endif

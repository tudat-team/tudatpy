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
#include "../../../shared_ptr_wasm.h"

#include <tudat/simulation/environment_setup/createGravityFieldVariations.h>

namespace tss = tudat::simulation_setup;
namespace tg = tudat::gravitation;

WASM_MODULE_PATH("dynamics_environment_setup_gravity_field_variation")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_environment_setup_gravity_field_variation) {
    using namespace emscripten;

    // BodyDeformationTypes enum
    enum_<tg::BodyDeformationTypes>(
        "dynamics_environment_setup_gravity_field_variation_BodyDeformationTypes")
        .value("basic_solid_body", tg::basic_solid_body)
        .value("iers_2010_tidal", tg::iers_2010)
        .value("tabulated_deformation", tg::tabulated_variation)
        .value("periodic_variation", tg::periodic_variation)
        .value("polynomial_variation", tg::polynomial_variation)
        .value("ocean_tide", tg::ocean_tide)
        .value("pole_tide", tg::pole_tide);

    // GravityFieldVariationSettings base class
    class_<tss::GravityFieldVariationSettings>(
        "dynamics_environment_setup_gravity_field_variation_GravityFieldVariationSettings")
        .smart_ptr<std::shared_ptr<tss::GravityFieldVariationSettings>>(
            "shared_ptr_GravityFieldVariationSettings");

    // BasicSolidBodyGravityFieldVariationSettings derived class
    class_<tss::BasicSolidBodyGravityFieldVariationSettings,
           base<tss::GravityFieldVariationSettings>>(
        "dynamics_environment_setup_gravity_field_variation_BasicSolidBodyGravityFieldVariationSettings")
        .smart_ptr<std::shared_ptr<tss::BasicSolidBodyGravityFieldVariationSettings>>(
            "shared_ptr_BasicSolidBodyGravityFieldVariationSettings");

    // ========================================================================
    // Factory functions
    // ========================================================================

    // Solid body tide with single degree Love number
    function("dynamics_environment_setup_gravity_field_variation_solid_body_tide",
        select_overload<std::shared_ptr<tss::GravityFieldVariationSettings>(
            const std::string, const double, const int)>(
                &tss::fixedSingleDegreeLoveNumberGravityFieldVariationSettings));

    // Solid body tide with complex Love number
    function("dynamics_environment_setup_gravity_field_variation_solid_body_tide_complex_k",
        select_overload<std::shared_ptr<tss::GravityFieldVariationSettings>(
            const std::string, const std::complex<double>, const int)>(
                &tss::fixedSingleDegreeLoveNumberGravityFieldVariationSettings));

    // Solid body tide with degree-variable Love numbers
    function("dynamics_environment_setup_gravity_field_variation_solid_body_tide_degree_variable_k",
        select_overload<std::shared_ptr<tss::GravityFieldVariationSettings>(
            const std::string, std::map<int, double>)>(
                &tss::fixedSingleDegreeLoveNumberGravityFieldVariationSettings));

    // Solid body tide with degree-variable complex Love numbers
    function("dynamics_environment_setup_gravity_field_variation_solid_body_tide_degree_variable_complex_k",
        select_overload<std::shared_ptr<tss::GravityFieldVariationSettings>(
            const std::string, std::map<int, std::complex<double>>)>(
                &tss::fixedSingleDegreeLoveNumberGravityFieldVariationSettings));

    // Mode-coupled solid body tide
    function("dynamics_environment_setup_gravity_field_variation_mode_coupled_solid_body_tide",
        &tss::modeCoupledSolidBodyGravityFieldVariationSettings);

    // Periodic gravity field variations
    function("dynamics_environment_setup_gravity_field_variation_periodic",
        &tss::periodicGravityFieldVariationsSettings);

    // Single-period periodic gravity field variations
    function("dynamics_environment_setup_gravity_field_variation_single_period_periodic",
        &tss::periodicGravityFieldVariationsSettingsSingleFrequency);

    // Polynomial gravity field variations
    function("dynamics_environment_setup_gravity_field_variation_polynomial",
        &tss::polynomialGravityFieldVariationsSettings);
}

#endif

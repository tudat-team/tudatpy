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

#include <tudat/simulation/environment_setup/createGravityField.h>

namespace tss = tudat::simulation_setup;
namespace tpc = tudat::physical_constants;

WASM_MODULE_PATH("dynamics_environment_setup_gravity_field")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_environment_setup_gravity_field) {
    using namespace emscripten;

    // GravityFieldType enum
    enum_<tss::GravityFieldType>("dynamics_environment_setup_gravity_field_GravityFieldType")
        .value("central", tss::central)
        .value("central_spice", tss::central_spice)
        .value("spherical_harmonic", tss::spherical_harmonic)
        .value("polyhedron", tss::polyhedron)
        .value("ring_gravity", tss::one_dimensional_ring);

    // SphericalHarmonicsModel enum (PredefinedSphericalHarmonicsModel in Python)
    enum_<tss::SphericalHarmonicsModel>("dynamics_environment_setup_gravity_field_SphericalHarmonicsModel")
        .value("customModel", tss::customModel)
        .value("egm96", tss::egm96)
        .value("ggm02c", tss::ggm02c)
        .value("ggm02s", tss::ggm02s)
        .value("goco05c", tss::goco05c)
        .value("glgm3150", tss::glgm3150)
        .value("lpe200", tss::lpe200)
        .value("gggrx1200", tss::gggrx1200)
        .value("jgmro120d", tss::jgmro120d)
        .value("jgmess160a", tss::jgmess160a)
        .value("shgj180u", tss::shgj180u);

    // GravityFieldSettings base class
    class_<tss::GravityFieldSettings>("dynamics_environment_setup_gravity_field_GravityFieldSettings")
        .smart_ptr<std::shared_ptr<tss::GravityFieldSettings>>("shared_ptr_GravityFieldSettings")
        .function("getGravityFieldType", &tss::GravityFieldSettings::getGravityFieldType);

    // CentralGravityFieldSettings
    class_<tss::CentralGravityFieldSettings, base<tss::GravityFieldSettings>>(
        "dynamics_environment_setup_gravity_field_CentralGravityFieldSettings")
        .smart_ptr<std::shared_ptr<tss::CentralGravityFieldSettings>>("shared_ptr_CentralGravityFieldSettings")
        .function("getGravitationalParameter", &tss::CentralGravityFieldSettings::getGravitationalParameter);

    // SphericalHarmonicsGravityFieldSettings
    class_<tss::SphericalHarmonicsGravityFieldSettings, base<tss::GravityFieldSettings>>(
        "dynamics_environment_setup_gravity_field_SphericalHarmonicsGravityFieldSettings")
        .smart_ptr<std::shared_ptr<tss::SphericalHarmonicsGravityFieldSettings>>(
            "shared_ptr_SphericalHarmonicsGravityFieldSettings")
        .function("getGravitationalParameter", &tss::SphericalHarmonicsGravityFieldSettings::getGravitationalParameter)
        .function("getReferenceRadius", &tss::SphericalHarmonicsGravityFieldSettings::getReferenceRadius)
        .function("getCosineCoefficients", &tss::SphericalHarmonicsGravityFieldSettings::getCosineCoefficients)
        .function("getSineCoefficients", &tss::SphericalHarmonicsGravityFieldSettings::getSineCoefficients);

    // FromFileSphericalHarmonicsGravityFieldSettings
    class_<tss::FromFileSphericalHarmonicsGravityFieldSettings, base<tss::SphericalHarmonicsGravityFieldSettings>>(
        "dynamics_environment_setup_gravity_field_FromFileSphericalHarmonicsGravityFieldSettings")
        .smart_ptr<std::shared_ptr<tss::FromFileSphericalHarmonicsGravityFieldSettings>>(
            "shared_ptr_FromFileSphericalHarmonicsGravityFieldSettings");

    // PolyhedronGravityFieldSettings
    class_<tss::PolyhedronGravityFieldSettings, base<tss::GravityFieldSettings>>(
        "dynamics_environment_setup_gravity_field_PolyhedronGravityFieldSettings")
        .smart_ptr<std::shared_ptr<tss::PolyhedronGravityFieldSettings>>(
            "shared_ptr_PolyhedronGravityFieldSettings")
        .function("getGravitationalParameter", &tss::PolyhedronGravityFieldSettings::getGravitationalParameter)
        .function("getDensity", &tss::PolyhedronGravityFieldSettings::getDensity);

    // Factory functions
    function("dynamics_environment_setup_gravity_field_central",
        &tss::centralGravitySettings);

    function("dynamics_environment_setup_gravity_field_central_spice",
        &tss::centralGravityFromSpiceSettings);

    function("dynamics_environment_setup_gravity_field_spherical_harmonic",
        select_overload<std::shared_ptr<tss::GravityFieldSettings>(
            const double, const double, const Eigen::MatrixXd, const Eigen::MatrixXd, const std::string&)>(
            &tss::sphericalHarmonicsGravitySettings));

    // Predefined spherical harmonics model
    function("dynamics_environment_setup_gravity_field_predefined_spherical_harmonic",
        optional_override([](tss::SphericalHarmonicsModel model, int maximumDegree) {
            return std::make_shared<tss::FromFileSphericalHarmonicsGravityFieldSettings>(
                model, maximumDegree);
        }));

    // Triaxial ellipsoid from density
    function("dynamics_environment_setup_gravity_field_sh_triaxial_ellipsoid_from_density",
        select_overload<std::shared_ptr<tss::SphericalHarmonicsGravityFieldSettings>(
            const double, const double, const double, const double,
            const int, const int, const std::string&, const double)>(
            &tss::createHomogeneousTriAxialEllipsoidGravitySettings));

    // Triaxial ellipsoid from gravitational parameter
    function("dynamics_environment_setup_gravity_field_sh_triaxial_ellipsoid_from_gravitational_parameter",
        select_overload<std::shared_ptr<tss::SphericalHarmonicsGravityFieldSettings>(
            const double, const double, const double,
            const int, const int, const std::string&, const double)>(
            &tss::createHomogeneousTriAxialEllipsoidGravitySettings));

    // Polyhedron from density
    function("dynamics_environment_setup_gravity_field_polyhedron_from_density",
        select_overload<std::shared_ptr<tss::GravityFieldSettings>(
            const double, const Eigen::MatrixXd, const Eigen::MatrixXi,
            const std::string&, const double)>(
            &tss::polyhedronGravitySettings));

    // Polyhedron from gravitational parameter
    function("dynamics_environment_setup_gravity_field_polyhedron_from_mu",
        select_overload<std::shared_ptr<tss::GravityFieldSettings>(
            const double, const Eigen::MatrixXd, const Eigen::MatrixXi,
            const std::string&, const double)>(
            &tss::polyhedronGravitySettingsFromMu));

    // Ring model
    function("dynamics_environment_setup_gravity_field_ring_model",
        &tss::ringGravitySettings);
}

#endif

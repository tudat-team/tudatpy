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
#include "../../../stl_wasm.h"
#include "../../../shared_ptr_wasm.h"

#include <tudat/simulation/environment_setup/createRadiationPressureInterface.h>

namespace tss = tudat::simulation_setup;

WASM_MODULE_PATH("dynamics_environment_setup_radiation_pressure")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_environment_setup_radiation_pressure) {
    using namespace emscripten;

    // RadiationPressureType enum
    enum_<tss::RadiationPressureType>("dynamics_environment_setup_radiation_pressure_RadiationPressureType")
        .value("cannon_ball_radiation_pressure_interface", tss::cannon_ball_radiation_pressure_interface);

    // RadiationPressureInterfaceSettings base class
    class_<tss::RadiationPressureInterfaceSettings>(
        "dynamics_environment_setup_radiation_pressure_RadiationPressureInterfaceSettings")
        .smart_ptr<std::shared_ptr<tss::RadiationPressureInterfaceSettings>>(
            "shared_ptr_RadiationPressureInterfaceSettings")
        .function("getRadiationPressureType", &tss::RadiationPressureInterfaceSettings::getRadiationPressureType)
        .function("getSourceBody", &tss::RadiationPressureInterfaceSettings::getSourceBody)
        .function("getOccultingBodies", &tss::RadiationPressureInterfaceSettings::getOccultingBodies);

    // CannonBallRadiationPressureInterfaceSettings
    class_<tss::CannonBallRadiationPressureInterfaceSettings, base<tss::RadiationPressureInterfaceSettings>>(
        "dynamics_environment_setup_radiation_pressure_CannonBallRadiationPressureInterfaceSettings")
        .smart_ptr<std::shared_ptr<tss::CannonBallRadiationPressureInterfaceSettings>>(
            "shared_ptr_CannonBallRadiationPressureInterfaceSettings")
        .function("getArea", &tss::CannonBallRadiationPressureInterfaceSettings::getArea)
        .function("getRadiationPressureCoefficient",
            &tss::CannonBallRadiationPressureInterfaceSettings::getRadiationPressureCoefficient);

    // Factory function for cannonball radiation pressure
    function("dynamics_environment_setup_radiation_pressure_cannonball",
        select_overload<std::shared_ptr<tss::RadiationPressureInterfaceSettings>(
            const std::string&, const double, const double, const std::vector<std::string>&)>(
            &tss::cannonBallRadiationPressureSettings));
}

#endif

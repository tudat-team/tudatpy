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

#include <tudat/astro/observation_models/linkTypeDefs.h>

namespace tom = tudat::observation_models;

WASM_MODULE_PATH("estimation_observable_models_setup_links")

EMSCRIPTEN_BINDINGS(tudatpy_estimation_observable_models_setup_links) {
    using namespace emscripten;

    // LinkEndType enum
    enum_<tom::LinkEndType>("estimation_observable_models_setup_links_LinkEndType")
        .value("unidentified_link_end", tom::unidentified_link_end)
        .value("transmitter", tom::transmitter)
        .value("reflector1", tom::reflector1)
        .value("reflector2", tom::reflector2)
        .value("reflector3", tom::reflector3)
        .value("reflector4", tom::reflector4)
        .value("receiver", tom::receiver)
        .value("observed_body", tom::observed_body)
        .value("transmitter2", tom::transmitter2)
        .value("observer", tom::observer);

    // LinkEndId struct
    class_<tom::LinkEndId>("estimation_observable_models_setup_links_LinkEndId")
        .constructor<>()
        .constructor<const std::string&>()
        .constructor<const std::string&, const std::string&>()
        .function("getBodyName", &tom::LinkEndId::getBodyName)
        .function("getStationName", &tom::LinkEndId::getStationName);

    // LinkDefinition struct
    class_<tom::LinkDefinition>("estimation_observable_models_setup_links_LinkDefinition")
        .constructor<>()
        .constructor<const std::map<tom::LinkEndType, tom::LinkEndId>&>()
        .function("getLinkEnds", &tom::LinkDefinition::getLinkEnds)
        .function("size", &tom::LinkDefinition::size);

    // Factory function for link end id
    function("estimation_observable_models_setup_links_link_end_id",
        select_overload<tom::LinkEndId(const std::string&, const std::string&)>(
            &tom::linkEndId));

    // Factory function for link definition
    function("estimation_observable_models_setup_links_link_definition",
        &tom::linkDefinition);
}

#endif

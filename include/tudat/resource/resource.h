/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Stub resource header for WebAssembly builds.
 *    In WASM builds, resource paths should be configured at runtime
 *    using Emscripten's virtual filesystem.
 *
 *    Users can mount data files to /tudat_data in the virtual filesystem.
 */

#ifndef TUDAT_RESOURCE_RESOURCE_H
#define TUDAT_RESOURCE_RESOURCE_H

#include <string>
#include <cstdlib>

// Default data path for WASM builds - can be configured via Emscripten's virtual FS
#ifndef TUDAT_WASM_DATA_PATH
#define TUDAT_WASM_DATA_PATH "/tudat_data"
#endif

namespace tudat
{
namespace paths
{

// Stub functions that provide path access for WASM builds
// These return default paths that can be populated via Emscripten's virtual filesystem

inline std::string get_resources_path()
{
    return TUDAT_WASM_DATA_PATH;
}

inline std::string get_ephemeris_path()
{
    return std::string(TUDAT_WASM_DATA_PATH) + "/ephemeris";
}

inline std::string get_earth_orientation_path()
{
    return std::string(TUDAT_WASM_DATA_PATH) + "/earth_orientation";
}

inline std::string get_quadrature_path()
{
    return std::string(TUDAT_WASM_DATA_PATH) + "/quadrature";
}

inline std::string get_spice_kernels_path()
{
    return std::string(TUDAT_WASM_DATA_PATH) + "/spice_kernels";
}

inline std::string get_atmosphere_tables_path()
{
    return std::string(TUDAT_WASM_DATA_PATH) + "/atmosphere_tables";
}

inline std::string get_gravity_models_path()
{
    return std::string(TUDAT_WASM_DATA_PATH) + "/gravity_models";
}

inline std::string get_space_weather_path()
{
    return std::string(TUDAT_WASM_DATA_PATH) + "/space_weather";
}

} // namespace paths
} // namespace tudat

#endif // TUDAT_RESOURCE_RESOURCE_H

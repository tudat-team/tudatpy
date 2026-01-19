# Tudatpy WASM Bindings CMake Configuration
# This creates Emscripten Embind bindings that mirror the PyBind11 API
#
# Build with:
#   emcmake cmake -B build-wasm -DCMAKE_BUILD_TYPE=Release ..
#   cmake --build build-wasm --target tudatpy_wasm

# Only build when Emscripten is detected
if(NOT TUDAT_BUILD_WASM)
    message(STATUS "Skipping WASM bindings (not an Emscripten build)")
    return()
endif()

message(STATUS "")
message(STATUS "*** CONFIGURING TUDATPY WASM BINDINGS ***")
message(STATUS "")

# ============================================================================
# Source Files
# ============================================================================

# Infrastructure headers (included by other files)
set(WASM_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/wasm_module.h
    ${CMAKE_CURRENT_SOURCE_DIR}/eigen_wasm.h
    ${CMAKE_CURRENT_SOURCE_DIR}/stl_wasm.h
    ${CMAKE_CURRENT_SOURCE_DIR}/shared_ptr_wasm.h
)

# Core binding source - main entry point
set(WASM_KERNEL_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/kernel_wasm.cpp
)

# Collect all module binding sources (will grow as modules are implemented)
set(WASM_MODULE_SOURCES "")

# Constants module
file(GLOB CONSTANTS_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/constants/*.cpp")
list(APPEND WASM_MODULE_SOURCES ${CONSTANTS_SOURCES})

# Math module
file(GLOB_RECURSE MATH_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/math/*.cpp")
list(APPEND WASM_MODULE_SOURCES ${MATH_SOURCES})

# Astro module
file(GLOB_RECURSE ASTRO_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/astro/*.cpp")
list(APPEND WASM_MODULE_SOURCES ${ASTRO_SOURCES})

# Dynamics module
file(GLOB_RECURSE DYNAMICS_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/dynamics/*.cpp")
list(APPEND WASM_MODULE_SOURCES ${DYNAMICS_SOURCES})

# Estimation module
file(GLOB_RECURSE ESTIMATION_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/estimation/*.cpp")
list(APPEND WASM_MODULE_SOURCES ${ESTIMATION_SOURCES})

# Interface module
file(GLOB_RECURSE INTERFACE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/interface/*.cpp")
list(APPEND WASM_MODULE_SOURCES ${INTERFACE_SOURCES})

# Data module
file(GLOB_RECURSE DATA_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/data/*.cpp")
list(APPEND WASM_MODULE_SOURCES ${DATA_SOURCES})

# Trajectory design module
file(GLOB_RECURSE TRAJECTORY_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/trajectory_design/*.cpp")
list(APPEND WASM_MODULE_SOURCES ${TRAJECTORY_SOURCES})

# Exceptions module
file(GLOB_RECURSE EXCEPTIONS_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/exceptions/*.cpp")
list(APPEND WASM_MODULE_SOURCES ${EXCEPTIONS_SOURCES})

# Combine all sources
set(ALL_WASM_SOURCES
    ${WASM_KERNEL_SOURCES}
    ${WASM_MODULE_SOURCES}
)

# ============================================================================
# Static Library Target (for linking into different WASM outputs)
# ============================================================================

add_library(tudatpy_wasm_bindings STATIC ${ALL_WASM_SOURCES})

target_include_directories(tudatpy_wasm_bindings
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_BINARY_DIR}/include
)

target_include_directories(tudatpy_wasm_bindings
    SYSTEM PUBLIC
        ${EIGEN3_INCLUDE_DIR}
        ${Boost_INCLUDE_DIRS}
)

# Add TUDAT_VERSION definition
target_compile_definitions(tudatpy_wasm_bindings PRIVATE
    TUDAT_VERSION="${PROJECT_VERSION}"
)

# Link against Tudat libraries (same dependencies as PyBind version)
# This MUST match Tudat_ESTIMATION_LIBRARIES for full 1:1 parity with native Tudat
target_link_libraries(tudatpy_wasm_bindings
    PUBLIC
        # Estimation setup (top-level)
        tudat_estimation_setup
        tudat_propagation_setup

        # Shape-based and trajectory design
        tudat_shape_based_methods
        tudat_low_thrust_trajectories

        # Environment setup
        tudat_environment_setup

        # Observation and estimation (CRITICAL for full estimation support)
        tudat_observation_models
        tudat_ground_stations
        tudat_acceleration_partials
        tudat_torque_partials
        tudat_observation_partials
        tudat_orbit_determination
        tudat_estimatable_parameters

        # Physics models
        tudat_aerodynamics
        tudat_system_models
        tudat_geometric_shapes
        tudat_relativity
        tudat_gravitation
        tudat_mission_segments
        tudat_electromagnetism
        tudat_propulsion
        tudat_ephemerides

        # Earth orientation (requires SOFA)
        tudat_earth_orientation

        # Core numerical methods
        tudat_numerical_integrators
        tudat_reference_frames
        tudat_statistics
        tudat_propagators

        # External interfaces
        tudat_sofa_interface
        tudat_spice_interface

        # Basic mathematics and utilities
        tudat_basic_astrodynamics
        tudat_numerical_quadrature
        tudat_interpolators
        tudat_root_finders
        tudat_basic_mathematics
        tudat_input_output
        tudat_basics
        tudat_data
)

# ============================================================================
# Main WASM Executable Target
# ============================================================================

add_executable(tudatpy_wasm ${WASM_KERNEL_SOURCES})

target_link_libraries(tudatpy_wasm
    PRIVATE
        -Wl,--whole-archive
        tudatpy_wasm_bindings
        -Wl,--no-whole-archive
)

# Emscripten-specific linker flags
set(WASM_LINK_FLAGS
    # Environment settings
    "-sENVIRONMENT=web,node"
    "-sMODULARIZE=1"
    "-sEXPORT_NAME='createTudatModule'"

    # Memory configuration
    "-sALLOW_MEMORY_GROWTH=1"
    "-sINITIAL_MEMORY=256MB"

    # Filesystem support (required for SPICE kernel loading)
    "-sFORCE_FILESYSTEM=1"

    # Function pointer emulation (required for f2c code in CSpice)
    "-sEMULATE_FUNCTION_POINTER_CASTS=1"

    # Embind support
    "--bind"

    # Exception handling
    "-sDISABLE_EXCEPTION_CATCHING=0"

    # Export additional runtime methods for JavaScript interop
    "-sEXPORTED_RUNTIME_METHODS=['cwrap','ccall','getValue','setValue','UTF8ToString','stringToUTF8','lengthBytesUTF8','FS']"
)

# Join flags with space
string(JOIN " " WASM_LINK_FLAGS_STRING ${WASM_LINK_FLAGS})

set_target_properties(tudatpy_wasm PROPERTIES
    SUFFIX ".js"
    LINK_FLAGS "${WASM_LINK_FLAGS_STRING}"
)

# ============================================================================
# Node.js Test Target (for automated testing)
# ============================================================================

add_executable(tudatpy_wasm_node ${WASM_KERNEL_SOURCES})

target_link_libraries(tudatpy_wasm_node
    PRIVATE
        tudatpy_wasm_bindings
)

set(WASM_NODE_LINK_FLAGS
    "-sENVIRONMENT=node"
    "-sMODULARIZE=1"
    "-sEXPORT_NAME='createTudatModule'"
    "-sALLOW_MEMORY_GROWTH=1"
    "-sINITIAL_MEMORY=256MB"
    "-sEMULATE_FUNCTION_POINTER_CASTS=1"
    "--bind"
    "-sDISABLE_EXCEPTION_CATCHING=0"
    "-sEXIT_RUNTIME=1"
    "-sNODERAWFS=1"
)

string(JOIN " " WASM_NODE_LINK_FLAGS_STRING ${WASM_NODE_LINK_FLAGS})

set_target_properties(tudatpy_wasm_node PROPERTIES
    SUFFIX ".js"
    LINK_FLAGS "${WASM_NODE_LINK_FLAGS_STRING}"
)

# ============================================================================
# Installation
# ============================================================================

install(TARGETS tudatpy_wasm tudatpy_wasm_node
    RUNTIME DESTINATION lib/wasm
)

install(FILES
    ${CMAKE_BINARY_DIR}/src/tudatpy_wasm/tudatpy_wasm.js
    ${CMAKE_BINARY_DIR}/src/tudatpy_wasm/tudatpy_wasm.wasm
    DESTINATION lib/wasm
    OPTIONAL
)

# ============================================================================
# NPM Package Target
# ============================================================================

set(NPM_PACKAGE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/npm)
set(NPM_DIST_DIR ${NPM_PACKAGE_DIR}/dist)

# Create npm package distribution directory
add_custom_target(npm_package_prepare
    COMMAND ${CMAKE_COMMAND} -E make_directory ${NPM_DIST_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/tudatpy_wasm.js
        ${NPM_DIST_DIR}/tudatpy_wasm.js
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/tudatpy_wasm.wasm
        ${NPM_DIST_DIR}/tudatpy_wasm.wasm
    DEPENDS tudatpy_wasm
    COMMENT "Preparing npm package distribution files"
)

# Generate ESM wrapper for modern JavaScript imports
add_custom_command(
    OUTPUT ${NPM_DIST_DIR}/tudatpy_wasm.mjs
    COMMAND ${CMAKE_COMMAND} -E echo "// ESM wrapper for Tudat WASM module" > ${NPM_DIST_DIR}/tudatpy_wasm.mjs
    COMMAND ${CMAKE_COMMAND} -E echo "import createTudatModuleFactory from './tudatpy_wasm.js';" >> ${NPM_DIST_DIR}/tudatpy_wasm.mjs
    COMMAND ${CMAKE_COMMAND} -E echo "export default createTudatModuleFactory;" >> ${NPM_DIST_DIR}/tudatpy_wasm.mjs
    COMMAND ${CMAKE_COMMAND} -E echo "export { createTudatModuleFactory as createTudatModule };" >> ${NPM_DIST_DIR}/tudatpy_wasm.mjs
    DEPENDS npm_package_prepare
    COMMENT "Generating ESM wrapper"
)

# Generate namespace builder for hierarchical API
configure_file(
    ${NPM_PACKAGE_DIR}/namespace-builder.js.in
    ${NPM_DIST_DIR}/tudatpy-namespace.js
    @ONLY
)

configure_file(
    ${NPM_PACKAGE_DIR}/namespace-builder.mjs.in
    ${NPM_DIST_DIR}/tudatpy-namespace.mjs
    @ONLY
)

# Main npm package target
add_custom_target(npm_package
    DEPENDS
        npm_package_prepare
        ${NPM_DIST_DIR}/tudatpy_wasm.mjs
    COMMENT "NPM package ready at ${NPM_PACKAGE_DIR}"
)

# NPM publish target (dry-run by default for safety)
add_custom_target(npm_publish_dry_run
    COMMAND npm pack --dry-run
    WORKING_DIRECTORY ${NPM_PACKAGE_DIR}
    DEPENDS npm_package
    COMMENT "Dry-run npm publish (use npm_publish for actual publish)"
)

# Actual npm publish (requires authentication)
add_custom_target(npm_publish
    COMMAND npm publish --access public
    WORKING_DIRECTORY ${NPM_PACKAGE_DIR}
    DEPENDS npm_package
    COMMENT "Publishing npm package to registry"
)

# ============================================================================
# Web Test Target
# ============================================================================

set(WEB_TEST_DIR ${CMAKE_SOURCE_DIR}/tests/wasm/web)

# Copy WASM files to web test directory for browser testing
add_custom_target(wasm_test_deploy
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/tudatpy_wasm.js
        ${WEB_TEST_DIR}/tudatpy_wasm.js
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/tudatpy_wasm.wasm
        ${WEB_TEST_DIR}/tudatpy_wasm.wasm
    DEPENDS tudatpy_wasm
    COMMENT "Deploying WASM files to web test directory"
)

# ============================================================================
# Documentation Target
# ============================================================================

find_program(NPM_EXECUTABLE npm)
find_program(NPX_EXECUTABLE npx)

if(NPM_EXECUTABLE AND NPX_EXECUTABLE)
    # Install documentation dependencies
    add_custom_target(docs_wasm_deps
        COMMAND ${NPM_EXECUTABLE} install
        WORKING_DIRECTORY ${NPM_PACKAGE_DIR}
        COMMENT "Installing documentation dependencies"
    )

    # Generate HTML documentation using TypeDoc
    add_custom_target(docs_wasm
        COMMAND ${NPX_EXECUTABLE} typedoc --options typedoc.json
        WORKING_DIRECTORY ${NPM_PACKAGE_DIR}
        DEPENDS docs_wasm_deps
        COMMENT "Generating WASM API documentation"
    )

    # Serve documentation locally for preview
    add_custom_target(docs_wasm_serve
        COMMAND ${NPX_EXECUTABLE} serve docs -l 3000
        WORKING_DIRECTORY ${NPM_PACKAGE_DIR}
        DEPENDS docs_wasm
        COMMENT "Serving documentation at http://localhost:3000"
    )

    message(STATUS "Documentation targets available: docs_wasm, docs_wasm_serve")
else()
    message(STATUS "npm/npx not found - documentation targets disabled")
endif()

# ============================================================================
# Status Messages
# ============================================================================

message(STATUS "WASM binding sources: ${ALL_WASM_SOURCES}")
message(STATUS "WASM output: tudatpy_wasm.js / tudatpy_wasm.wasm")
message(STATUS "NPM package directory: ${NPM_PACKAGE_DIR}")
message(STATUS "")
message(STATUS "Available WASM targets:")
message(STATUS "  tudatpy_wasm        - Build web+node WASM module")
message(STATUS "  tudatpy_wasm_node   - Build node-only WASM module")
message(STATUS "  wasm_test_deploy    - Deploy WASM to tests/wasm/web/ for browser testing")
message(STATUS "  npm_package         - Prepare npm package for publishing")
message(STATUS "  npm_publish_dry_run - Test npm publish without uploading")
message(STATUS "  npm_publish         - Publish to npm registry")
message(STATUS "  docs_wasm           - Generate HTML documentation")
message(STATUS "  docs_wasm_serve     - Serve documentation locally")
message(STATUS "")

# WASM-specific tests that don't require Boost.Test
# These tests can be run with Node.js after building with Emscripten
#
# This test suite includes comprehensive propagation tests that work without
# external SPICE kernels by using analytical ephemerides and self-contained
# dynamical systems (CR3BP, two-body, mass propagation, etc.)

# Source files for the WASM test suite
set(WASM_TEST_SOURCES
    src/main.cpp
    src/testBasicAstro.cpp
    src/testPropagation.cpp
    src/testSpice.cpp
    src/testGravitation.cpp
    src/testAerodynamics.cpp
    src/testMissionSegments.cpp
    src/testIntegrators.cpp
    src/testElectromagnetism.cpp
    src/testEphemerides.cpp
    src/testEarthOrientation.cpp
    src/testExamples.cpp
    src/testEstimation.cpp
    src/testEdgeCases.cpp
)

# Emscripten bindings (only for web build with --bind flag)
set(WASM_BINDINGS_SOURCES
    src/wasmBindings.cpp
)

add_executable(tudat_wasm_test ${WASM_TEST_SOURCES})

target_include_directories(tudat_wasm_test
    PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(tudat_wasm_test
    SYSTEM PRIVATE
    ${EIGEN3_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
)

# Full propagation library dependencies for comprehensive WASM testing
# These libraries enable orbit propagation, environment setup, and dynamics simulation
# without requiring external SPICE kernels (using analytical ephemerides instead)
target_link_libraries(tudat_wasm_test
    PRIVATE
    # Simulation and propagation setup
    tudat_propagation_setup
    tudat_environment_setup
    tudat_estimation_setup

    # Core propagators and dynamics
    tudat_propagators
    tudat_numerical_integrators

    # Astrodynamics and physics
    tudat_basic_astrodynamics
    tudat_gravitation
    tudat_ephemerides
    tudat_aerodynamics
    tudat_propulsion
    tudat_electromagnetism
    tudat_relativity
    tudat_earth_orientation

    # Supporting libraries
    tudat_reference_frames
    tudat_basic_mathematics
    tudat_interpolators
    tudat_root_finders
    tudat_statistics
    tudat_numerical_quadrature

    # System and geometry
    tudat_system_models
    tudat_geometric_shapes
    tudat_ground_stations

    # Mission analysis
    tudat_mission_segments
    tudat_shape_based_methods
    tudat_low_thrust_trajectories

    # I/O and utilities
    tudat_input_output
    tudat_basics
    tudat_data

    # External interfaces (SOFA for high-precision calculations, SPICE for ephemeris)
    tudat_sofa_interface
    tudat_spice_interface
    ${Sofa_LIBRARIES}
    ${CSpice_LIBRARIES}

    # Atmosphere models
    ${NRLMSISE00_LIBRARIES}

    # CALCEPH for direct binary SPK reading (WASM only)
    ${CALCEPH_LIBRARIES}
)

set_target_properties(tudat_wasm_test
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/tests/wasm"
)

# For Emscripten, we need to configure for Node.js execution
if(EMSCRIPTEN)
    # Embed data files into the WASM binary
    # The data directory contains EOP files needed for earth orientation calculations
    set(WASM_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")

    # -sEMULATE_FUNCTION_POINTER_CASTS=1 is required because CSPICE's f2c-generated code
    # has function signature mismatches (e.g., s_cat, s_copy return void but are called
    # expecting int). WASM is strict about this and will crash without emulation.
    set_target_properties(tudat_wasm_test PROPERTIES
        SUFFIX ".js"
        LINK_FLAGS "-sENVIRONMENT=node -sEXIT_RUNTIME=1 -sNODERAWFS=0 -sMODULARIZE=0 -sINVOKE_RUN=1 -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=256MB -sEMULATE_FUNCTION_POINTER_CASTS=1 --embed-file ${WASM_DATA_DIR}@/tudat_data"
    )

    # Register test with CTest for running via Node.js
    # The test output is tudat_wasm_test.js in the build directory
    find_program(NODE_EXECUTABLE NAMES node nodejs)
    if(NODE_EXECUTABLE)
        add_test(
            NAME tudat_wasm_test
            COMMAND ${NODE_EXECUTABLE} "${PROJECT_BINARY_DIR}/tests/wasm/tudat_wasm_test.js"
            WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/tests/wasm"
        )
        set_tests_properties(tudat_wasm_test PROPERTIES
            LABELS "wasm"
        )
    else()
        message(WARNING "Node.js not found - WASM tests cannot be run via CTest. Install Node.js or run manually.")
    endif()
endif()

# ============================================================================
# Web-based test runner with UI
# ============================================================================

if(EMSCRIPTEN)
    # Web version of the test suite (runs in browser)
    # Includes Emscripten bindings for exposing C++ functions to JavaScript
    add_executable(tudat_wasm_web ${WASM_TEST_SOURCES} ${WASM_BINDINGS_SOURCES})

    target_include_directories(tudat_wasm_web
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_BINARY_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src
    )

    target_include_directories(tudat_wasm_web
        SYSTEM PRIVATE
        ${EIGEN3_INCLUDE_DIR}
        ${Boost_INCLUDE_DIRS}
    )

    target_link_libraries(tudat_wasm_web
        PRIVATE
        tudat_propagation_setup
        tudat_environment_setup
        tudat_estimation_setup
        tudat_propagators
        tudat_numerical_integrators
        tudat_basic_astrodynamics
        tudat_gravitation
        tudat_ephemerides
        tudat_aerodynamics
        tudat_propulsion
        tudat_electromagnetism
        tudat_relativity
        tudat_earth_orientation
        tudat_reference_frames
        tudat_basic_mathematics
        tudat_interpolators
        tudat_root_finders
        tudat_statistics
        tudat_numerical_quadrature
        tudat_system_models
        tudat_geometric_shapes
        tudat_ground_stations
        tudat_mission_segments
        tudat_shape_based_methods
        tudat_low_thrust_trajectories
        tudat_input_output
        tudat_basics
        tudat_data
        tudat_sofa_interface
        tudat_spice_interface
        ${Sofa_LIBRARIES}
        ${CSpice_LIBRARIES}
        ${NRLMSISE00_LIBRARIES}
        # CALCEPH for direct binary SPK reading (WASM only)
        ${CALCEPH_LIBRARIES}
    )

    # Web-specific output directory
    set(WASM_WEB_DIR "${PROJECT_BINARY_DIR}/tests/wasm/web")
    file(MAKE_DIRECTORY ${WASM_WEB_DIR})

    set_target_properties(tudat_wasm_web PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${WASM_WEB_DIR}"
        SUFFIX ".js"
        OUTPUT_NAME "tudat_wasm_test"
        # Web-compatible flags:
        # - ENVIRONMENT=web for browser execution
        # - MODULARIZE=0 so it runs immediately
        # - INVOKE_RUN=0 so we can control when main() runs
        # - EXPORTED_RUNTIME_METHODS for JS interop
        # - lembind for Emscripten bindings (exposes C++ functions to JS)
        LINK_FLAGS "-sENVIRONMENT=web,node,worker -sEXIT_RUNTIME=0 -sNODERAWFS=0 -sMODULARIZE=1 -sEXPORT_NAME='createTudatModule' -sINVOKE_RUN=0 -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=256MB -sEMULATE_FUNCTION_POINTER_CASTS=1 -sFORCE_FILESYSTEM=1 -sEXPORTED_RUNTIME_METHODS=['callMain','FS','cwrap','ccall'] -sEXPORTED_FUNCTIONS=['_main'] --bind --embed-file ${WASM_DATA_DIR}@/tudat_data"
    )

    # Copy web UI files to build directory (with proper dependency tracking)
    # Using OUTPUT-based custom command so changes to source files trigger rebuild
    add_custom_command(
        OUTPUT "${WASM_WEB_DIR}/index.html" "${WASM_WEB_DIR}/app.js" "${WASM_WEB_DIR}/testWorker.js"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/web/index.html"
            "${WASM_WEB_DIR}/index.html"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/web/app.js"
            "${WASM_WEB_DIR}/app.js"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/web/testWorker.js"
            "${WASM_WEB_DIR}/testWorker.js"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/web/index.html"
            "${CMAKE_CURRENT_SOURCE_DIR}/web/app.js"
            "${CMAKE_CURRENT_SOURCE_DIR}/web/testWorker.js"
        COMMENT "Copying web UI files to build directory"
    )

    # Copy imagery files (Blue Marble and Black Marble for day/night Earth)
    file(MAKE_DIRECTORY "${WASM_WEB_DIR}/imagery")
    add_custom_command(
        OUTPUT "${WASM_WEB_DIR}/imagery/world.topo.bathy.200407.3x5400x2700.jpg"
               "${WASM_WEB_DIR}/imagery/BlackMarble_2016_3km.jpeg"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/web/imagery/world.topo.bathy.200407.3x5400x2700.jpg"
            "${WASM_WEB_DIR}/imagery/world.topo.bathy.200407.3x5400x2700.jpg"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/web/imagery/BlackMarble_2016_3km.jpeg"
            "${WASM_WEB_DIR}/imagery/BlackMarble_2016_3km.jpeg"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/web/imagery/world.topo.bathy.200407.3x5400x2700.jpg"
            "${CMAKE_CURRENT_SOURCE_DIR}/web/imagery/BlackMarble_2016_3km.jpeg"
        COMMENT "Copying Earth imagery files to build directory"
    )

    # Copy CesiumJS library (local copy for offline use)
    add_custom_command(
        OUTPUT "${WASM_WEB_DIR}/cesium/Cesium.js"
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${CMAKE_CURRENT_SOURCE_DIR}/web/cesium"
            "${WASM_WEB_DIR}/cesium"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/web/cesium/Cesium.js"
        COMMENT "Copying CesiumJS library to build directory"
    )

    # Copy visualizations directory (includes all visualization modules and examples)
    add_custom_command(
        OUTPUT "${WASM_WEB_DIR}/visualizations/index.js"
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${CMAKE_CURRENT_SOURCE_DIR}/web/visualizations"
            "${WASM_WEB_DIR}/visualizations"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/web/visualizations/index.js"
        COMMENT "Copying visualizations directory to build directory"
    )

    # Copy API documentation directory
    add_custom_command(
        OUTPUT "${WASM_WEB_DIR}/api/index.html"
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${CMAKE_CURRENT_SOURCE_DIR}/web/api"
            "${WASM_WEB_DIR}/api"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/web/api/index.html"
        COMMENT "Copying API documentation to build directory"
    )

    # Copy precomputed ephemeris data (from pre-converted SPK files)
    file(MAKE_DIRECTORY "${WASM_WEB_DIR}/data/ephemeris")
    add_custom_command(
        OUTPUT "${WASM_WEB_DIR}/data/ephemeris/earth_sun_j2000.json"
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${CMAKE_CURRENT_SOURCE_DIR}/web/data/ephemeris"
            "${WASM_WEB_DIR}/data/ephemeris"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/web/data/ephemeris/earth_sun_j2000.json"
        COMMENT "Copying precomputed ephemeris data to build directory"
    )

    # Create a target for the web files that tudat_wasm_web depends on
    add_custom_target(tudat_wasm_web_files
        DEPENDS "${WASM_WEB_DIR}/index.html" "${WASM_WEB_DIR}/app.js" "${WASM_WEB_DIR}/testWorker.js"
                "${WASM_WEB_DIR}/imagery/world.topo.bathy.200407.3x5400x2700.jpg"
                "${WASM_WEB_DIR}/imagery/BlackMarble_2016_3km.jpeg"
                "${WASM_WEB_DIR}/cesium/Cesium.js"
                "${WASM_WEB_DIR}/visualizations/index.js"
                "${WASM_WEB_DIR}/api/index.html"
                "${WASM_WEB_DIR}/data/ephemeris/earth_sun_j2000.json"
    )
    add_dependencies(tudat_wasm_web tudat_wasm_web_files)

    # ============================================================================
    # Web server target: starts a local server and opens the browser
    # ============================================================================

    # Find Python for the web server
    find_package(Python3 COMPONENTS Interpreter)

    if(Python3_FOUND)
        # Create a script to start the server
        file(WRITE "${WASM_WEB_DIR}/start_server.py" [=[
#!/usr/bin/env python3
"""
Simple HTTP server for Tudat WASM test suite.
Serves files with correct MIME types for WASM.
"""
import http.server
import socketserver
import webbrowser
import os
import sys

PORT = 8832  # TUD on phone keypad (T=8, U=8, D=3) + 2 for "to"
DIRECTORY = os.path.dirname(os.path.abspath(__file__))

class ReusableTCPServer(socketserver.TCPServer):
    """TCP server that allows address reuse for clean restarts."""
    allow_reuse_address = True

class WasmHandler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=DIRECTORY, **kwargs)

    def end_headers(self):
        # Required headers for WASM with SharedArrayBuffer support
        # Using 'credentialless' instead of 'require-corp' to allow CDN resources
        self.send_header('Cross-Origin-Opener-Policy', 'same-origin')
        self.send_header('Cross-Origin-Embedder-Policy', 'credentialless')
        super().end_headers()

    def guess_type(self, path):
        if path.endswith('.wasm'):
            return 'application/wasm'
        if path.endswith('.js'):
            return 'application/javascript'
        return super().guess_type(path)

def main():
    os.chdir(DIRECTORY)

    with ReusableTCPServer(("", PORT), WasmHandler) as httpd:
        url = f"http://localhost:{PORT}"
        print(f"\n{'='*60}")
        print(f"  Tudat WASM Test Suite")
        print(f"{'='*60}")
        print(f"\n  Server running at: {url}")
        print(f"  Press Ctrl+C to stop\n")
        print(f"{'='*60}\n")

        # Open browser
        webbrowser.open(url)

        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\n\nServer stopped.")
            sys.exit(0)

if __name__ == "__main__":
    main()
]=])

        # Add custom target to run the web server
        add_custom_target(tudat_wasm_serve
            COMMAND ${Python3_EXECUTABLE} "${WASM_WEB_DIR}/start_server.py"
            DEPENDS tudat_wasm_web
            WORKING_DIRECTORY "${WASM_WEB_DIR}"
            COMMENT "Starting Tudat WASM test server at http://localhost:8080"
            USES_TERMINAL
        )

        message(STATUS "Web server target available: ninja tudat_wasm_serve")

        # Combined target: build, test, and serve
        # This is a convenience target that runs the Node.js tests first, then starts the web server
        if(NODE_EXECUTABLE)
            add_custom_target(wasm
                COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --target tudat_wasm_test
                COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --target tudat_wasm_web
                COMMAND ${CMAKE_COMMAND} -E echo "Running WASM tests via Node.js..."
                COMMAND ${NODE_EXECUTABLE} "${PROJECT_BINARY_DIR}/tests/wasm/tudat_wasm_test.js"
                COMMAND ${CMAKE_COMMAND} -E echo "Tests passed! Starting web server..."
                COMMAND ${Python3_EXECUTABLE} "${WASM_WEB_DIR}/start_server.py"
                WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
                COMMENT "Building, testing, and serving Tudat WASM"
                USES_TERMINAL
            )
            message(STATUS "Combined WASM target available: ninja wasm")
        else()
            # If Node.js not available, just build and serve without testing
            add_custom_target(wasm
                COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --target tudat_wasm_web
                COMMAND ${CMAKE_COMMAND} -E echo "Starting web server (tests skipped - Node.js not found)..."
                COMMAND ${Python3_EXECUTABLE} "${WASM_WEB_DIR}/start_server.py"
                WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
                COMMENT "Building and serving Tudat WASM (tests skipped)"
                USES_TERMINAL
            )
            message(STATUS "Combined WASM target available: ninja wasm (tests disabled - Node.js not found)")
        endif()
    else()
        message(WARNING "Python3 not found - web server target not available")
    endif()

endif()

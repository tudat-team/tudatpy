# Tudat WASM build, test, and serve script
# Usage: cmake -P wasm.cmake
#        cmake -DFULL_BUILD=1 -P wasm.cmake   (includes full API + tests)
#
# This script automatically downloads Emscripten SDK if needed, then
# configures, builds, tests, and serves the WASM version of Tudat.

cmake_minimum_required(VERSION 3.20)

# Option for full build (can be set via -DFULL_BUILD=1)
option(FULL_BUILD "Build full tudatpy_wasm API and run tests" OFF)

# Get the directory containing this script (the tudat root)
get_filename_component(TUDAT_ROOT "${CMAKE_CURRENT_LIST_DIR}" ABSOLUTE)

set(WASM_BUILD_DIR "${TUDAT_ROOT}/build-wasm")
set(EMSDK_DIR "${TUDAT_ROOT}/.emsdk")
set(EMSDK_INSTALL_DIR "${EMSDK_DIR}/emsdk")
set(EMSDK_VERSION "3.1.51")
set(CUSTOM_TOOLCHAIN "${TUDAT_ROOT}/cmake_modules/toolchain-emscripten.cmake")

# Determine emsdk executable path
if(WIN32)
    set(EMSDK_EXE "${EMSDK_INSTALL_DIR}/emsdk.bat")
else()
    set(EMSDK_EXE "${EMSDK_INSTALL_DIR}/emsdk")
endif()

message(STATUS "")
message(STATUS "============================================================")
message(STATUS "  Tudat WASM Build, Test, and Serve")
message(STATUS "============================================================")
message(STATUS "")

# ============================================================
# Step 1: Download and install Emscripten SDK if needed
# ============================================================
if(NOT EXISTS "${EMSDK_INSTALL_DIR}")
    message(STATUS "Emscripten SDK not found. Downloading...")
    file(MAKE_DIRECTORY "${EMSDK_DIR}")

    execute_process(
        COMMAND git clone --depth 1 https://github.com/emscripten-core/emsdk.git "${EMSDK_INSTALL_DIR}"
        RESULT_VARIABLE GIT_RESULT
    )
    if(NOT GIT_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to clone Emscripten SDK")
    endif()
    message(STATUS "Emscripten SDK cloned successfully")
endif()

# Check if the requested version is installed
set(EMSDK_VERSION_FILE "${EMSDK_DIR}/.installed_version")
set(NEED_INSTALL TRUE)

if(EXISTS "${EMSDK_VERSION_FILE}")
    file(READ "${EMSDK_VERSION_FILE}" INSTALLED_VERSION)
    string(STRIP "${INSTALLED_VERSION}" INSTALLED_VERSION)
    if("${INSTALLED_VERSION}" STREQUAL "${EMSDK_VERSION}")
        set(NEED_INSTALL FALSE)
        message(STATUS "Emscripten ${EMSDK_VERSION} is already installed")
    else()
        message(STATUS "Emscripten version change: ${INSTALLED_VERSION} -> ${EMSDK_VERSION}")
    endif()
endif()

if(NEED_INSTALL)
    message(STATUS "Installing Emscripten ${EMSDK_VERSION}...")

    execute_process(
        COMMAND "${EMSDK_EXE}" install ${EMSDK_VERSION}
        WORKING_DIRECTORY "${EMSDK_INSTALL_DIR}"
        RESULT_VARIABLE INSTALL_RESULT
    )
    if(NOT INSTALL_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to install Emscripten ${EMSDK_VERSION}")
    endif()

    execute_process(
        COMMAND "${EMSDK_EXE}" activate ${EMSDK_VERSION}
        WORKING_DIRECTORY "${EMSDK_INSTALL_DIR}"
        RESULT_VARIABLE ACTIVATE_RESULT
    )
    if(NOT ACTIVATE_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to activate Emscripten ${EMSDK_VERSION}")
    endif()

    file(WRITE "${EMSDK_VERSION_FILE}" "${EMSDK_VERSION}")
    message(STATUS "Emscripten ${EMSDK_VERSION} installed and activated")
endif()

# ============================================================
# Step 2: Configure the build
# ============================================================
file(MAKE_DIRECTORY "${WASM_BUILD_DIR}")

if(EXISTS "${WASM_BUILD_DIR}/CMakeCache.txt")
    message(STATUS "Build already configured, skipping configure step...")
    message(STATUS "(Delete build-wasm/CMakeCache.txt to force reconfigure)")
else()
    message(STATUS "Configuring WASM build...")
    execute_process(
        COMMAND ${CMAKE_COMMAND} -G Ninja
            -DCMAKE_BUILD_TYPE=Release
            -DCMAKE_TOOLCHAIN_FILE=${CUSTOM_TOOLCHAIN}
            "${TUDAT_ROOT}"
        WORKING_DIRECTORY "${WASM_BUILD_DIR}"
        RESULT_VARIABLE CONFIG_RESULT
    )
    if(NOT CONFIG_RESULT EQUAL 0)
        message(FATAL_ERROR "CMake configuration failed")
    endif()
endif()

# ============================================================
# Step 3: Build WASM modules
# ============================================================

# Build the visualization test module (tudat_wasm_web) - includes propagation visualizations
message(STATUS "")
message(STATUS "Building tudat_wasm_web (visualization demo)...")
execute_process(
    COMMAND ${CMAKE_COMMAND} --build . --target tudat_wasm_web -j 8
    WORKING_DIRECTORY "${WASM_BUILD_DIR}"
    RESULT_VARIABLE BUILD_RESULT
)
if(NOT BUILD_RESULT EQUAL 0)
    message(WARNING "tudat_wasm_web build failed (may not be configured)")
endif()

# Skip full Embind API build by default (slow) - set FULL_BUILD=1 to include it
if(DEFINED ENV{FULL_BUILD} OR FULL_BUILD)
    message(STATUS "")
    message(STATUS "Building tudatpy_wasm (full Embind API)...")
    execute_process(
        COMMAND ${CMAKE_COMMAND} --build . --target tudatpy_wasm -j 8
        WORKING_DIRECTORY "${WASM_BUILD_DIR}"
        RESULT_VARIABLE BUILD_RESULT
    )
    if(NOT BUILD_RESULT EQUAL 0)
        message(WARNING "tudatpy_wasm build failed")
    endif()
else()
    message(STATUS "")
    message(STATUS "Skipping tudatpy_wasm build (set FULL_BUILD=1 to include)")
endif()

# ============================================================
# Step 4: Deploy to test directories
# ============================================================
message(STATUS "")
message(STATUS "Deploying WASM files...")

set(EMBIND_OUTPUT_DIR "${WASM_BUILD_DIR}/src/tudatpy_wasm")
set(LEGACY_OUTPUT_DIR "${WASM_BUILD_DIR}/tests/wasm/web")
set(WEB_TEST_DIR "${TUDAT_ROOT}/tests/wasm/web")
set(NPM_DIST_DIR "${TUDAT_ROOT}/src/tudatpy_wasm/npm/dist")

# Deploy visualization module to web test directory
file(MAKE_DIRECTORY "${WEB_TEST_DIR}")
if(EXISTS "${LEGACY_OUTPUT_DIR}/tudat_wasm_test.js")
    file(COPY "${LEGACY_OUTPUT_DIR}/tudat_wasm_test.js" DESTINATION "${WEB_TEST_DIR}")
    file(COPY "${LEGACY_OUTPUT_DIR}/tudat_wasm_test.wasm" DESTINATION "${WEB_TEST_DIR}")
    message(STATUS "  -> Deployed tudat_wasm_test to tests/wasm/web/")
endif()

# Deploy visualization module to docs directory (for GitHub Pages)
set(DOCS_DIR "${TUDAT_ROOT}/docs")
if(EXISTS "${LEGACY_OUTPUT_DIR}/tudat_wasm_test.js")
    file(COPY "${LEGACY_OUTPUT_DIR}/tudat_wasm_test.js" DESTINATION "${DOCS_DIR}")
    file(COPY "${LEGACY_OUTPUT_DIR}/tudat_wasm_test.wasm" DESTINATION "${DOCS_DIR}")
    message(STATUS "  -> Deployed tudat_wasm_test to docs/")
endif()

# Deploy full Embind API module to web test directory (only if built)
if(EXISTS "${EMBIND_OUTPUT_DIR}/tudatpy_wasm.js")
    file(COPY "${EMBIND_OUTPUT_DIR}/tudatpy_wasm.js" DESTINATION "${WEB_TEST_DIR}")
    file(COPY "${EMBIND_OUTPUT_DIR}/tudatpy_wasm.wasm" DESTINATION "${WEB_TEST_DIR}")
    message(STATUS "  -> Deployed tudatpy_wasm to tests/wasm/web/")
endif()

# Deploy Embind module to npm package dist directory
file(MAKE_DIRECTORY "${NPM_DIST_DIR}")
if(EXISTS "${EMBIND_OUTPUT_DIR}/tudatpy_wasm.js")
    file(COPY "${EMBIND_OUTPUT_DIR}/tudatpy_wasm.js" DESTINATION "${NPM_DIST_DIR}")
    file(COPY "${EMBIND_OUTPUT_DIR}/tudatpy_wasm.wasm" DESTINATION "${NPM_DIST_DIR}")
    message(STATUS "  -> Deployed tudatpy_wasm to src/tudatpy_wasm/npm/dist/")
endif()

# ============================================================
# Step 5: Run Node.js tests (only if tudatpy_wasm was built)
# ============================================================
if(EXISTS "${EMBIND_OUTPUT_DIR}/tudatpy_wasm.js")
    message(STATUS "")
    message(STATUS "Running WASM tests via Node.js...")
    find_program(NODE_EXECUTABLE NAMES node nodejs)
    if(NODE_EXECUTABLE)
        # Run the API test
        if(EXISTS "${TUDAT_ROOT}/tests/wasm/test_wasm_api.js")
            execute_process(
                COMMAND ${NODE_EXECUTABLE} "${TUDAT_ROOT}/tests/wasm/test_wasm_api.js"
                WORKING_DIRECTORY "${TUDAT_ROOT}/tests/wasm"
                RESULT_VARIABLE TEST_RESULT
            )
            if(TEST_RESULT EQUAL 0)
                message(STATUS "API tests passed!")
            endif()
        endif()

        # Run the integration test
        if(EXISTS "${TUDAT_ROOT}/tests/wasm/test_wasm_integration.js")
            execute_process(
                COMMAND ${NODE_EXECUTABLE} "${TUDAT_ROOT}/tests/wasm/test_wasm_integration.js"
                WORKING_DIRECTORY "${TUDAT_ROOT}/tests/wasm"
                RESULT_VARIABLE TEST_RESULT
            )
            if(TEST_RESULT EQUAL 0)
                message(STATUS "Integration tests passed!")
            endif()
        endif()
    endif()
endif()

# ============================================================
# Step 6: Start web server
# ============================================================
message(STATUS "")
message(STATUS "============================================================")
message(STATUS "  Starting web server at http://localhost:8832")
message(STATUS "  Press Ctrl+C to stop")
message(STATUS "============================================================")
message(STATUS "")

find_package(Python3 COMPONENTS Interpreter REQUIRED)

# Create a simple server script in the web test directory
file(WRITE "${WEB_TEST_DIR}/start_server.py" [=[
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
        print(f"  Tudat WASM Embind Test Suite")
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

execute_process(
    COMMAND ${Python3_EXECUTABLE} "${WEB_TEST_DIR}/start_server.py"
    WORKING_DIRECTORY "${WEB_TEST_DIR}"
)

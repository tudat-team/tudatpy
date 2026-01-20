#    Copyright (c) 2010-2024, Delft University of Technology
#    All rights reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#
#    Emscripten Toolchain File for WebAssembly builds
#
#    This toolchain automatically downloads and installs the Emscripten SDK
#    if not already present. The SDK is installed to .emsdk/ in the project root.
#
#    Usage:
#      cmake -B build-wasm -DCMAKE_TOOLCHAIN_FILE=cmake_modules/toolchain-emscripten.cmake
#      cmake --build build-wasm
#
#    To update Emscripten:
#      cmake --build build-wasm --target update-emscripten
#

# Emscripten version to use
set(EMSDK_VERSION "3.1.51" CACHE STRING "Emscripten SDK version to install")

# Find git before any environment changes - needed for FetchContent
find_program(GIT_EXECUTABLE git)
if(GIT_EXECUTABLE)
    set(CMAKE_GIT_EXECUTABLE "${GIT_EXECUTABLE}" CACHE FILEPATH "Git executable" FORCE)
endif()

# Get the project root directory (where this toolchain file lives is cmake_modules/)
get_filename_component(TUDAT_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
set(EMSDK_DIR "${TUDAT_ROOT}/.emsdk")
set(EMSDK_INSTALL_DIR "${EMSDK_DIR}/emsdk")

# Function to download and setup Emscripten SDK
function(setup_emscripten_sdk)
    if(NOT EXISTS "${EMSDK_INSTALL_DIR}")
        message(STATUS "")
        message(STATUS "==============================================")
        message(STATUS "  Emscripten SDK not found. Installing...")
        message(STATUS "==============================================")
        message(STATUS "")

        # Create directory
        file(MAKE_DIRECTORY "${EMSDK_DIR}")

        # Clone emsdk repository
        message(STATUS "Cloning Emscripten SDK...")
        execute_process(
            COMMAND git clone --depth 1 https://github.com/emscripten-core/emsdk.git "${EMSDK_INSTALL_DIR}"
            RESULT_VARIABLE GIT_RESULT
            OUTPUT_VARIABLE GIT_OUTPUT
            ERROR_VARIABLE GIT_ERROR
        )
        if(NOT GIT_RESULT EQUAL 0)
            message(FATAL_ERROR "Failed to clone Emscripten SDK: ${GIT_ERROR}")
        endif()
    endif()

    # Check if the requested version is installed and active
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

        # Determine emsdk executable
        if(WIN32)
            set(EMSDK_EXE "${EMSDK_INSTALL_DIR}/emsdk.bat")
        else()
            set(EMSDK_EXE "${EMSDK_INSTALL_DIR}/emsdk")
        endif()

        # Install the specified version
        execute_process(
            COMMAND "${EMSDK_EXE}" install ${EMSDK_VERSION}
            WORKING_DIRECTORY "${EMSDK_INSTALL_DIR}"
            RESULT_VARIABLE INSTALL_RESULT
            OUTPUT_VARIABLE INSTALL_OUTPUT
            ERROR_VARIABLE INSTALL_ERROR
        )
        if(NOT INSTALL_RESULT EQUAL 0)
            message(FATAL_ERROR "Failed to install Emscripten ${EMSDK_VERSION}: ${INSTALL_ERROR}\n${INSTALL_OUTPUT}")
        endif()

        # Activate the version
        execute_process(
            COMMAND "${EMSDK_EXE}" activate ${EMSDK_VERSION}
            WORKING_DIRECTORY "${EMSDK_INSTALL_DIR}"
            RESULT_VARIABLE ACTIVATE_RESULT
            OUTPUT_VARIABLE ACTIVATE_OUTPUT
            ERROR_VARIABLE ACTIVATE_ERROR
        )
        if(NOT ACTIVATE_RESULT EQUAL 0)
            message(FATAL_ERROR "Failed to activate Emscripten ${EMSDK_VERSION}: ${ACTIVATE_ERROR}\n${ACTIVATE_OUTPUT}")
        endif()

        # Save installed version
        file(WRITE "${EMSDK_VERSION_FILE}" "${EMSDK_VERSION}")

        message(STATUS "Emscripten ${EMSDK_VERSION} installed and activated")
    endif()
endfunction()

# Run setup
setup_emscripten_sdk()

# Set up paths to Emscripten tools
set(EMSCRIPTEN_ROOT "${EMSDK_INSTALL_DIR}/upstream/emscripten")

# System identification
set(CMAKE_SYSTEM_NAME Emscripten)
set(CMAKE_SYSTEM_VERSION 1)
set(EMSCRIPTEN TRUE CACHE BOOL "Emscripten build")

# Set compilers
set(CMAKE_C_COMPILER "${EMSCRIPTEN_ROOT}/emcc")
set(CMAKE_CXX_COMPILER "${EMSCRIPTEN_ROOT}/em++")
set(CMAKE_AR "${EMSCRIPTEN_ROOT}/emar" CACHE FILEPATH "Emscripten ar")
set(CMAKE_RANLIB "${EMSCRIPTEN_ROOT}/emranlib" CACHE FILEPATH "Emscripten ranlib")

# Mark compilers as working
set(CMAKE_C_COMPILER_WORKS TRUE)
set(CMAKE_CXX_COMPILER_WORKS TRUE)

# Emscripten uses Clang under the hood
set(CMAKE_C_COMPILER_ID Clang)
set(CMAKE_CXX_COMPILER_ID Clang)

# Library suffixes
set(CMAKE_STATIC_LIBRARY_SUFFIX ".a")
set(CMAKE_EXECUTABLE_SUFFIX ".js")

# Position independent code (always on for Emscripten)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Compile flags (no -s flags here - those are linker settings)
set(CMAKE_C_FLAGS_INIT "")
set(CMAKE_CXX_FLAGS_INIT "")

# Build type specific flags
set(CMAKE_C_FLAGS_DEBUG_INIT "-g -O0")
set(CMAKE_C_FLAGS_RELEASE_INIT "-O3 -DNDEBUG")
set(CMAKE_C_FLAGS_RELWITHDEBINFO_INIT "-O2 -g -DNDEBUG")
set(CMAKE_C_FLAGS_MINSIZEREL_INIT "-Os -DNDEBUG")

set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -g -DNDEBUG")
set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT "-Os -DNDEBUG")

# Linker flags for executables/modules
# -s WASM=1: Output WebAssembly (default in recent Emscripten)
# -s DISABLE_EXCEPTION_CATCHING=0: Enable C++ exceptions (required for Tudat)
# -s ALLOW_MEMORY_GROWTH=1: Allow heap to grow dynamically
# Note: MODULARIZE and EXPORT_ES6 are set per-target since they vary by use case
set(CMAKE_EXE_LINKER_FLAGS_INIT "-s WASM=1 -s DISABLE_EXCEPTION_CATCHING=0 -s ALLOW_MEMORY_GROWTH=1")

# Force appropriate Tudat build options for WASM
# Tests require Boost.Test which doesn't work well in WASM
set(TUDAT_BUILD_TESTS OFF CACHE BOOL "Disabled for WASM build" FORCE)

# Tutorials typically produce executables, disable for library build
set(TUDAT_BUILD_TUDAT_TUTORIALS OFF CACHE BOOL "Disabled for WASM build" FORCE)

# Static libraries are the norm for WASM
set(TUDAT_BUILD_STATIC_LIBRARY ON CACHE BOOL "Static library for WASM" FORCE)

# SOFA interface is required for WASM builds due to unconditional includes in source
# It will be automatically fetched and built with Emscripten
set(TUDAT_BUILD_WITH_SOFA_INTERFACE ON CACHE BOOL "Required for WASM build" FORCE)

# Extended precision not well supported in WASM
set(TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS OFF CACHE BOOL "Disabled for WASM" FORCE)

# Find root configuration for cross-compilation
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)

# Use Emscripten's Boost port
# This provides Boost headers compiled for WebAssembly
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -sUSE_BOOST_HEADERS=1")

# Pre-fetch Emscripten's Boost port so CMake can find it
# The embuilder tool downloads and builds Emscripten ports
message(STATUS "Fetching Emscripten Boost port...")
execute_process(
    COMMAND "${EMSCRIPTEN_ROOT}/embuilder" build boost_headers
    WORKING_DIRECTORY "${EMSCRIPTEN_ROOT}"
    RESULT_VARIABLE BOOST_FETCH_RESULT
    OUTPUT_VARIABLE BOOST_FETCH_OUTPUT
    ERROR_VARIABLE BOOST_FETCH_ERROR
)
if(NOT BOOST_FETCH_RESULT EQUAL 0)
    message(WARNING "Failed to fetch Boost port: ${BOOST_FETCH_ERROR}")
endif()


# Set Boost variables to help CMake find the Emscripten-provided Boost
# The port installs headers to the sysroot cache
set(BOOST_ROOT "${EMSDK_INSTALL_DIR}/upstream/emscripten/cache/sysroot")
set(Boost_INCLUDE_DIR "${EMSDK_INSTALL_DIR}/upstream/emscripten/cache/sysroot/include")
set(Boost_NO_SYSTEM_PATHS ON)
set(Boost_NO_BOOST_CMAKE ON)

# Disable CMake package registry to avoid finding stale Eigen3 entries
# from other projects in /private/tmp or other locations that may not exist
set(CMAKE_FIND_PACKAGE_NO_PACKAGE_REGISTRY ON CACHE BOOL "Disable package registry for clean WASM builds" FORCE)

cmake_minimum_required(VERSION 2.8.3)

if (POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

#CMP0048 is set to NEW
project(tudatpy VERSION 0.3.1 LANGUAGES CXX C)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/external/cmake")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}-Wl,--no-undefined")

message(STATUS "System name: ${CMAKE_SYSTEM_NAME}")
message(STATUS "tudatpy version: ${tudatpy_VERSION}")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# `PyLong_Type' can not be used when making a shared object; recompile with -fPIC
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# Load UserSettings.txt
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    set(BUILD_STYLE "standalone")
    include("${CMAKE_CURRENT_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
else ()
    set(BUILD_STYLE "part of ${CMAKE_PROJECT_NAME}")
    include("${CMAKE_CURRENT_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
    include("${CMAKE_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
    STRING(REGEX REPLACE ${CMAKE_SOURCE_DIR} "" RELATIVE_PROJECT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
    set(RELATIVE_PROJECT_PATH "${RELATIVE_PROJECT_PATH}" CACHE STRING "Relative path wrt to project for function")
    # message(STATUS "Relative path (wrt to project): ${RELATIVE_PROJECT_PATH}")
endif ()


# Set CMake build-type. If it not supplied by the user (either directly as an argument of through
# the "UserSettings.txt" file, the default built type is "Release".
if ((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
    set(CMAKE_BUILD_TYPE Release)
elseif (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_BUILD_TYPE Debug)
endif ()

message(STATUS "<< ${PROJECT_NAME} (${CMAKE_BUILD_TYPE} - ${BUILD_STYLE}) >>")

# Add local module path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMakeModules")
message(STATUS "CMake Module path(s): ${CMAKE_MODULE_PATH}")

# Set compiler based on preferences (e.g. USE_CLANG) and system.
include(compiler)

# Define the directory with the source code.
set(SRCROOT "${CMAKE_CURRENT_SOURCE_DIR}/tudatpy")

# Define the code root directory.
set(CODEROOT "${CMAKE_CURRENT_SOURCE_DIR}")

# Set testing options based on platform.
enable_testing()

# Set lib and bin directories where static libraries and unit tests are built.
if (NOT LIBROOT)
    set(LIBROOT "${CODEROOT}/lib")
endif ()
if (NOT BINROOT)
    set(BINROOT "${CODEROOT}/bin")
endif ()

# Set the global macros for setting up targets.
macro(setup_tudatpy_library_target target_name CUSTOM_OUTPUT_PATH)
    set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIBROOT}")
    set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIBROOT}")
endmacro(setup_tudatpy_library_target)

# Include the top-level directories.
include_directories(AFTER
        "${CODEROOT}"
        )

# Find Eigen3 library on local system.
find_package(Eigen3 REQUIRED)

# Include Eigen3 directories.
# Set CMake flag to suppress Eigen warnings (platform-dependent solution).
if (NOT APPLE OR APPLE_INCLUDE_FORCE)
    include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")
else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
endif ()

# Configure Boost libraries.
if (NOT Boost_USE_STATIC_LIBS)
    set(Boost_USE_STATIC_LIBS ON)
endif ()
if (NOT Boost_USE_MULTITHREADED)
    set(Boost_USE_MULTITHREADED ON)
endif ()
if (NOT Boost_USE_STATIC_RUNTIME)
    set(Boost_USE_STATIC_RUNTIME ON)
endif ()

# Find Boost libraries on local system.
find_package(Boost 1.55.0
        COMPONENTS thread date_time system unit_test_framework filesystem regex python3 REQUIRED)

# Include Boost directories.
# Set CMake flag to suppress Boost warnings (platform-dependent solution).
if (NOT APPLE OR APPLE_INCLUDE_FORCE)
    include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")
else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif ()

# Find Tudat library on local system.
find_package(Tudat 2.0 REQUIRED)

# Include Tudat directories.
# Set CMake flag to suppress Tudat warnings (platform-dependent solution).
if (NOT APPLE OR APPLE_INCLUDE_FORCE)
    include_directories(SYSTEM AFTER "${TUDAT_INCLUDE_DIR}")
else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${TUDAT_INCLUDE_DIR}\"")
endif ()

# Find CSPICE library on local system.
find_package(Spice)

# Include CSpice directories.
if (NOT APPLE OR APPLE_INCLUDE_FORCE)
    include_directories(SYSTEM AFTER "${SPICE_INCLUDE_DIR}")
else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${SPICE_INCLUDE_DIR}\"")
endif ()

if (NOT USE_NRLMSISE00)
    message(STATUS "NRLMSISE-00 disabled!")
    add_definitions(-DUSE_NRLMSISE00=0)
else ()
    message(STATUS "NRLMSISE-00 enabled!")
    add_definitions(-DUSE_NRLMSISE00=1)
    # Find USE_NRLMSISE00 library on local system.
    find_package(NRLMSISE00)

    # Include NRLMSISE00 directories.
    if (NOT APPLE OR APPLE_INCLUDE_FORCE)
        include_directories(SYSTEM AFTER "${NRLMSISE00_INCLUDE_DIR}")
    else ()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${NRLMSISE00_INCLUDE_DIR}\"")
    endif ()
endif ()

# Set compiler based on preferences (e.g. USE_CLANG) and system.
include(tudatLinkLibraries)

set(POSITION_INDEPENDENT_CODE True)

enable_testing()

# python.
include(YACMAPythonSetup)

# python version check.
if (${PYTHON_VERSION_MAJOR} LESS 3 OR (${PYTHON_VERSION_MAJOR} EQUAL 3 AND ${PYTHON_VERSION_MINOR} LESS 4))
    message(FATAL_ERROR "Minimum supported python version is 3.4.")
endif ()

# pybind11.
find_package(pybind11 REQUIRED)

# Configure the sphinx config file.
#configure_file("${CMAKE_CURRENT_SOURCE_DIR}/doc/conf.py.in" "${CMAKE_CURRENT_SOURCE_DIR}/doc/conf.py" @ONLY)

# Add the module directory.
add_subdirectory(tudatpy)

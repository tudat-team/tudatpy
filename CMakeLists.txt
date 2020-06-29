#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.

# Specify minimum cmake version required to compile Tudat.
cmake_minimum_required(VERSION 3.3)

# Trick for recompiling changes for BUILD-TREE.
if (EXISTS "${CMAKE_BINARY_DIR}/CMakeCache.txt")
    file(REMOVE "${CMAKE_BINARY_DIR}/CMakeCache.txt")
endif ()

# Set cmake policy.
cmake_policy(SET CMP0069 NEW)

# Get version from file.
file(READ "version" Tudat_VERSION)

# Specific project name.
project(Tudat VERSION ${Tudat_VERSION} LANGUAGES C CXX)

# Get lower string name of PROJECT_NAME for directory configuration.
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWER)
string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPER)
set(CMAKE_POSITION_INDEPENDENT_CODE True)

# Add local module path.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/yolo")

# General system information.
message(STATUS "System name: ${CMAKE_SYSTEM_NAME}")
enable_testing()
set(Tudat_DEFINITIONS "")

# +============================================================================
# BUILD OPTIONS
#  Offer the user the choice of defining the build variation.
# +============================================================================
# Build option: enable the test suite.
option(TUDAT_BUILD_TESTS "Build the test suite." ON)

# Build option: include default data suite.
option(TUDAT_FETCH_DATA_SUITE "Downloads default data suite during build." OFF)

# Build with propagation tests.
option(TUDAT_BUILD_WITH_PROPAGATION_TESTS "Build tudat with propagation tests. (>30 s propagations - Total test time > 10 minutes.)" OFF)

# Build with estimation tools.
option(TUDAT_BUILD_WITH_ESTIMATION_TOOLS "Build tudat with estimation tools." OFF)

# Build option: enable tutorials.
option(TUDAT_BUILD_TUDAT_TUTORIALS "Build tutorials." ON)

# Build option: enable tutorials.
option(TUDAT_BUILD_PAGMO_TUTORIALS "Build pagmo2 tutorials." OFF)

# Build option: enable tutorials.
option(TUDAT_BUILD_TBB_SUPPORT "Build tbb support for Pagmo." OFF)

# Build static library instead of dynamic.
option(TUDAT_BUILD_STATIC_LIBRARY "Build tudat as a static library, instead of dynamic." ON)

# Build static library instead of dynamic.
option(TUDAT_BUILD_WITH_FILTERS "Build Tudat with preliminary Kalman Filter functionality." OFF)

# Build sofa interface.
option(TUDAT_BUILD_WITH_SOFA_INTERFACE "Build Tudat with sofa interface." ON)

# Build json interface.
option(TUDAT_BUILD_WITH_JSON_INTERFACE "Build Tudat with json interface." OFF)

# Build spice interface.
option(TUDAT_BUILD_WITH_SPICE_INTERFACE "Build Tudat with spice interface." ON)

# Build extended precision propagation tools.
option(TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS "Build tudat with extended precision propagation tools." OFF)

# +============================================================================
# INSTALL TREE CONFIGURATION (Project name independent)
#  Offer the user the choice of overriding the installation directories.
# +============================================================================
# TODO: Consider if this can work with <tudat/paths.hpp> going forward.
set(INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib/" CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin/" CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include/" CACHE PATH "Installation directory for header files")
set(INSTALL_DATA_DIR "${CMAKE_INSTALL_PREFIX}/data/" CACHE PATH "Installation directory for data files")
set(INSTALL_MAN_DIR "${CMAKE_INSTALL_PREFIX}/" CACHE PATH "Installation directory for man documentation")
#set(INSTALL_TESTS_DIR "${CMAKE_INSTALL_PREFIX}/tests/" CACHE PATH "Installation directory for tests (default=OFF)")
set(INSTALL_CMAKE_DIR "${INSTALL_LIB_DIR}/cmake/${PROJECT_NAME_LOWER}" CACHE PATH "Installation directory for cmake config files")

# Make relative paths absolute (needed later on)
foreach (p LIB BIN INCLUDE DATA CMAKE)
    set(var INSTALL_${p}_DIR)
    set(RELATIVE_INSTALL_${p}_DIR ${INSTALL_${p}_DIR})
    if (NOT IS_ABSOLUTE "${${var}}")
        set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
    endif ()
endforeach ()

# Find relative paths for cmake-config file. (for find_package)
file(RELATIVE_PATH rel_include_dir "${INSTALL_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")
file(RELATIVE_PATH rel_lib_dir "${INSTALL_CMAKE_DIR}" "${INSTALL_LIB_DIR}")
file(RELATIVE_PATH rel_data_dir "${INSTALL_CMAKE_DIR}" "${INSTALL_DATA_DIR}")
file(RELATIVE_PATH rel_bin_dir "${INSTALL_CMAKE_DIR}" "${INSTALL_BIN_DIR}")
file(RELATIVE_PATH rel_man_dir "${INSTALL_CMAKE_DIR}" "${INSTALL_MAN_DIR}")

# Set relative paths for config.cmake.
foreach (p include lib data bin man)
    string(TOUPPER ${p} P)
    set(RELATIVE_INSTALL_${P}_DIR ${rel_${p}_dir})
endforeach ()

# +============================================================================
# DEPENDENCIES SPECIFIC TO BUILD OPTIONS
# +============================================================================
if (NOT TUDAT_BUILD_STATIC_LIBRARY)
    # https://stackoverflow.com/questions/30398238/cmake-rpath-not-working-could-not-find-shared-object-file
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif ()

# Enforce TBB if pagmo tutorials are to be built.
if (TUDAT_BUILD_PAGMO_TUTORIALS)
    if (TUDAT_BUILD_TBB_SUPPORT)
    else ()
        message("Forcing TBB support for Pagmo tutorial build.")
        set(TUDAT_BUILD_TBB_SUPPORT ON)
    endif ()
endif ()

# Set cmake build-type. If it not supplied by the user (either directly as an argument of through
# the "UserSettings.txt" file, the default built type is "Release".
if ((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
    set(CMAKE_BUILD_TYPE Release)
elseif (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_BUILD_TYPE Debug)
endif ()

message(STATUS "<< ${PROJECT_NAME} (${CMAKE_BUILD_TYPE} - ${BUILD_STYLE}) >>")

# Set compiler based on preferences (e.g. USE_CLANG) and system.
include(compiler)

#if (TUDAT_CONDA_BUILD)
message(STATUS "Building conda-build build of tudat.")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -DTUDAT_CONDA_BUILD=${TUDAT_CONDA_BUILD} -D_GLIBCXX_USE_CXX11_ABI=0")
#endif ()

# Get tudat resource manager dependency.
find_package(TudatResources REQUIRED)
include_directories(TudatResources_INCLUDE_DIRS)

# Find Eigen3 library on local system and include directories.
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")
if (APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
endif ()

# Adds Boost_UNIT_TEST_FRAMEWORK_LIBRARY if required.
# Configure Boost libraries.
if (TUDAT_BUILD_TESTS)
    set(_TUDAT_FIND_BOOST_UNIT_TEST_FRAMEWORK ON)
endif ()
include(TudatFindBoost)
if (APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif ()

# CSpice dependency if in build settings.
if (TUDAT_BUILD_WITH_SPICE_INTERFACE)
    add_definitions(-DTUDAT_BUILD_WITH_SPICE_INTERFACE=1)
    if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
        find_package(CSpice REQUIRED 1.0.0)
    elseif (NOT ${TUDAT_BUNDLE_BUILD})
        find_package(CSpice REQUIRED 1.0.0)
    else ()
        set(CSpice_INCLUDE_DIRS ${CSpice_BUILD_INCLUDE_DIRS})
    endif ()
else ()
    add_definitions(-DTUDAT_BUILD_WITH_SPICE_INTERFACE=0)
endif ()

# Sofa dependency if in build settings.
if (TUDAT_BUILD_WITH_SOFA_INTERFACE)
    #
    add_definitions(-DTUDAT_BUILD_WITH_SOFA_INTERFACE=1)
    if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
        find_package(Sofa REQUIRED 1.0.0)
    elseif (NOT ${TUDAT_BUNDLE_BUILD})
        find_package(Sofa REQUIRED 1.0.0)
    else ()
        set(Sofa_INCLUDE_DIRS ${Sofa_BUILD_INCLUDE_DIRS})
    endif ()
else ()
    add_definitions(-DTUDAT_BUILD_WITH_SOFA_INTERFACE=0)
endif ()

# JSON
if (TUDAT_BUILD_WITH_JSON_INTERFACE)
    find_package(nlohmann_json REQUIRED 3.7.3)
endif ()

# NRLMSISE-00
if (NOT BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS)
    message(STATUS "Extended precision propagation disabled!")
    add_definitions(-DBUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=0)
else ()
    message(STATUS "Extended precision propagation enabled!")
    add_definitions(-DBUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=1)
endif ()

include(YOLOProjectAddTestCase)
include(YOLOProjectAddLibrary)
include(YOLOProjectAddExecutable)
include(YOLOProjectAddExternalData)

# Set compiler based on preferences (e.g. USE_CLANG) and system.
include(TudatLinkLibraries)

# Add to parent scope if project is added via add_subdirectory().
if (NOT ${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    set(Tudat_PROPAGATION_LIBRARIES ${Tudat_PROPAGATION_LIBRARIES} PARENT_SCOPE)
    set(Tudat_ESTIMATION_LIBRARIES ${Tudat_ESTIMATION_LIBRARIES} PARENT_SCOPE)
endif ()

add_subdirectory(src/astro)
add_subdirectory(src/basics)
add_subdirectory(src/interface)
add_subdirectory(src/math)
add_subdirectory(src/simulation)
add_subdirectory(src/io)
add_subdirectory(examples)

if (TUDAT_BUILD_TESTS)
    add_subdirectory(tests)
endif ()

# Cleanup YOLO global project variables.
#include(YOLOProjectCleanup)

file(RELATIVE_PATH TUDAT_DATA_DIR_RELATIVE_TO_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" "${INSTALL_DATA_DIR}/tudat")

#configure_file( # Really dirty trick, and should be removed.
#        "${PROJECT_SOURCE_DIR}/paths.hpp.in"
#        "${PROJECT_BINARY_DIR}/include/tudat/paths.hpp" @ONLY
#        )

configure_file(
        "${PROJECT_SOURCE_DIR}/config.hpp.in"
        "${PROJECT_BINARY_DIR}/include/tudat/config.hpp" @ONLY
)

# Installation of the header files.
install(DIRECTORY include/ DESTINATION "${INSTALL_INCLUDE_DIR}")
install(FILES
        "${PROJECT_BINARY_DIR}/include/tudat/config.hpp"
        DESTINATION ${INSTALL_INCLUDE_DIR}/tudat)

# Add all targets to the build-tree export set
export(TARGETS ${${PROJECT_NAME}_LIBRARIES} NAMESPACE ${PROJECT_NAME}:: FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}_export.cmake")

# Installation of the export modules.
install(EXPORT ${PROJECT_NAME_LOWER}_export NAMESPACE ${PROJECT_NAME}:: DESTINATION "${RELATIVE_INSTALL_CMAKE_DIR}")

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE ${PROJECT_NAME})

# +============================================================================
# PROJECT CONFIG
#  For external use of project.
# +============================================================================
include(CMakePackageConfigHelpers)

configure_package_config_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME_LOWER}-config.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}-config.cmake"
        INSTALL_DESTINATION "${INSTALL_CMAKE_DIR}")

configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME_LOWER}-config.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}-config.cmake" @ONLY
)

# Configure file path for tudat data loading.
#set(Tudat_DATA_DIRS ${INSTALL_CMAKE_DIR}/${RELATIVE_INSTALL_DATA_DIR})


file(
        COPY "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/TudatFindBoost.cmake"
        DESTINATION "${CMAKE_CURRENT_BINARY_DIR}"
)

install(
        FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}-config.cmake"
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/TudatFindBoost.cmake"
        DESTINATION "${INSTALL_CMAKE_DIR}"
)

# Take care of versioning.

# NOTE: SameMinorVersion available only
# since CMake 3.11.
if (${CMAKE_VERSION} VERSION_LESS "3.11.0")
    write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}-config-version.cmake" VERSION ${${PROJECT_NAME}_VERSION}
            COMPATIBILITY SameMajorVersion)
else ()
    write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}-config-version.cmake" VERSION ${${PROJECT_NAME}_VERSION}
            COMPATIBILITY SameMinorVersion)
endif ()
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}-config-version.cmake" DESTINATION "${INSTALL_CMAKE_DIR}")

# ========
message(STATUS "TUDAT_DATA_DIR_RELATIVE_TO_INSTALL_PREFIX: ${TUDAT_DATA_DIR_RELATIVE_TO_INSTALL_PREFIX}")

# +============================================================================
# CLEAN UP (Project name independent)
#  Cleanup in case project is not top level.
# +============================================================================
unset(INSTALL_LIB_DIR)
unset(INSTALL_BIN_DIR)
unset(INSTALL_INCLUDE_DIR)
unset(INSTALL_DATA_DIR)
unset(INSTALL_MAN_DIR)
unset(INSTALL_TESTS_DIR)
unset(INSTALL_CMAKE_DIR)

unset(RELATIVE_INSTALL_LIB_DIR)
unset(RELATIVE_INSTALL_BIN_DIR)
unset(RELATIVE_INSTALL_INCLUDE_DIR)
unset(RELATIVE_INSTALL_DATA_DIR)
unset(RELATIVE_INSTALL_MAN_DIR)
unset(RELATIVE_INSTALL_TESTS_DIR)
unset(RELATIVE_INSTALL_CMAKE_DIR)

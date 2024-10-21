#    Copyright (c) 2010-2022, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.

# TODO: https://dev.azure.com/tudat-team/feedstock-builds/_build/results?buildId=26&view=logs&j=4b3757f9-99c5-5891-6c0c-93d7030cdf7d&t=f047df50-a233-5ae9-bf8f-c7dd8e5684ad -> https://stackoverflow.com/questions/36567072/why-do-i-get-ld-warning-direct-access-in-main-to-global-weak-symbol-in-this

# Specify minimum cmake version required to compile Tudat.
cmake_minimum_required(VERSION 3.3)

# Trick for recompiling changes for BUILD-TREE.
if (EXISTS "${CMAKE_BINARY_DIR}/CMakeCache.txt")
    file(REMOVE "${CMAKE_BINARY_DIR}/CMakeCache.txt")
endif ()

# Set cmake policy.
cmake_policy(SET CMP0069 NEW)

# Get version from file.
file(READ "version" Tudat_VERSION_LABEL)
string(REGEX MATCH "^([0-9]+)\\.([0-9]+)\\.([0-9]+)" Tudat_VERSION "${Tudat_VERSION_LABEL}")

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

option(TUDAT_DOWNLOAD_AND_BUILD_BOOST "Downloads and builds boost" OFF)

# Build option: include default data suite.
#option(TUDAT_FETCH_DATA_SUITE "Downloads default data suite during build." OFF)

# Build with propagation tests.
option(TUDAT_BUILD_WITH_PROPAGATION_TESTS "Build tudat with propagation tests. (>30 s propagations - Total test time > 10 minutes.)" OFF)

# Build with estimation tools.
option(TUDAT_BUILD_WITH_ESTIMATION_TOOLS "Build tudat with estimation tools." ON)

# Build option: enable tutorials.
option(TUDAT_BUILD_TUDAT_TUTORIALS "Build tutorials." ON)

# Build static library instead of dynamic.
option(TUDAT_BUILD_STATIC_LIBRARY "Build tudat as a static library, instead of dynamic." ON)

# Build static library instead of dynamic.
option(TUDAT_BUILD_WITH_FILTERS "Build Tudat with preliminary Kalman Filter functionality." OFF)
 
# Build sofa interface.
option(TUDAT_BUILD_WITH_SOFA_INTERFACE "Build Tudat with sofa interface." ON)

# Build json interface.
option(TUDAT_BUILD_WITH_JSON_INTERFACE "Build Tudat with json interface." OFF)

option(TUDAT_BUILD_WITH_FFTW3 "Build Tudat with FFTW3." OFF)

# Build pagmo-dependent code
option(TUDAT_BUILD_WITH_PAGMO "Build Tudat with pagmo." OFF)
if(CMAKE_CXX_SIMULATE_ID MATCHES "MSVC")
    # Build extended precision propagation tools.
    option(TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS "Build tudat with extended precision propagation tools." OFF)
else()
    option(TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS "Build tudat with extended precision propagation tools." OFF)
endif()

# Build as part of a GitHub Actions workflow
# Option enables the use of ccache by MSVC
# see https://github.com/ccache/ccache/wiki/MS-Visual-Studio
option(TUDAT_BUILD_GITHUB_ACTIONS "Build as part of GitHub Actions workflow." OFF)
if (TUDAT_BUILD_GITHUB_ACTIONS AND MSVC)
    find_program(ccache_exe ccache)
    if(ccache_exe)
        message(STATUS "CONFIGURING MSVC FOR CCACHE")
        file(COPY_FILE ${ccache_exe} ${CMAKE_BINARY_DIR}/cl.exe ONLY_IF_DIFFERENT)

        # By default Visual Studio generators will use /Zi which is not compatible
        # with ccache, so tell Visual Studio to use /Z7 instead.
        message(STATUS "Setting MSVC debug information format to 'Embedded'")
        set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<$<CONFIG:Debug,RelWithDebInfo>:Embedded>")

        set(CMAKE_VS_GLOBALS
        "CLToolExe=cl.exe"
        "CLToolPath=${CMAKE_BINARY_DIR}"
        "TrackFileAccess=false"
        "UseMultiToolTask=true"
        "DebugInformationFormat=OldStyle"
        )
    endif()
endif()

message(STATUS "******************** BUILD CONFIGURATION ********************")
message(STATUS "TUDAT_BUILD_TESTS                                     ${TUDAT_BUILD_TESTS}")
message(STATUS "TUDAT_BUILD_WITH_PROPAGATION_TESTS                    ${TUDAT_BUILD_WITH_PROPAGATION_TESTS}")
message(STATUS "TUDAT_BUILD_WITH_ESTIMATION_TOOLS                     ${TUDAT_BUILD_WITH_ESTIMATION_TOOLS}")
message(STATUS "TUDAT_BUILD_TUDAT_TUTORIALS                           ${TUDAT_BUILD_TUDAT_TUTORIALS}")
message(STATUS "TUDAT_BUILD_STATIC_LIBRARY                            ${TUDAT_BUILD_STATIC_LIBRARY}")
message(STATUS "TUDAT_BUILD_WITH_FILTERS                              ${TUDAT_BUILD_WITH_FILTERS}")
message(STATUS "TUDAT_BUILD_WITH_SOFA_INTERFACE                       ${TUDAT_BUILD_WITH_SOFA_INTERFACE}")
message(STATUS "TUDAT_BUILD_WITH_FFTW3                                ${TUDAT_BUILD_WITH_FFTW3}")
message(STATUS "TUDAT_BUILD_WITH_JSON_INTERFACE                       ${TUDAT_BUILD_WITH_JSON_INTERFACE}")
message(STATUS "TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS ${TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS}")
message(STATUS "TUDAT_DOWNLOAD_AND_BUILD_BOOST                        ${TUDAT_DOWNLOAD_AND_BUILD_BOOST}")

set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -DTUDAT_BUILD_WITH_FILTERS=${TUDAT_BUILD_WITH_FILTERS}")
set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -DTUDAT_BUILD_WITH_SOFA_INTERFACE=${TUDAT_BUILD_WITH_SOFA_INTERFACE}")
set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -DTUDAT_BUILD_WITH_FFTW3=${TUDAT_BUILD_WITH_FFTW3}")
set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -DTUDAT_BUILD_WITH_JSON_INTERFACE=${TUDAT_BUILD_WITH_JSON_INTERFACE}")
set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -DTUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=${TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS}")
# +============================================================================
# INSTALL TREE CONFIGURATION (Project name independent)
#  Offer the user the choice of overriding the installation directories.
# +============================================================================
# TODO: Consider if this can work with <tudat/paths.hpp> going forward.

set(INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib")
set(INSTALL_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin")
set(INSTALL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include")
set(INSTALL_DATA_DIR "${CMAKE_INSTALL_PREFIX}/data")
set(INSTALL_MAN_DIR "${CMAKE_INSTALL_PREFIX}")
#set(INSTALL_TESTS_DIR "${CMAKE_INSTALL_PREFIX}/tests")
set(INSTALL_CMAKE_DIR "${INSTALL_LIB_DIR}/cmake/${PROJECT_NAME_LOWER}")

message("")
message(STATUS "INSTALLATION PREFIX: ${CMAKE_INSTALL_PREFIX}")
message("")


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

# Set required boost components.
set(_TUDAT_REQUIRED_BOOST_LIBS filesystem system regex date_time thread
        chrono atomic)

if (TUDAT_BUILD_TESTS)
    # Add the unit test framework, if needed.
    list(APPEND _TUDAT_REQUIRED_BOOST_LIBS unit_test_framework)
endif ()

# Download boost if option is used.
if (TUDAT_DOWNLOAD_AND_BUILD_BOOST OR (TUDAT_BOOST_VERSION AND TUDAT_BOOST_SHA1))
    set(TUDAT_DOWNLOAD_AND_BUILD_BOOST ON)
    # Set the required boost components (for tudat)
    set(BoostComponents ${_TUDAT_REQUIRED_BOOST_LIBS})
    # Pick a boost version
    if (TUDAT_BOOST_VERSION AND TUDAT_BOOST_SHA1)
        set(BoostVersion ${TUDAT_BOOST_VERSION})
        set(BoostSHA1 ${TUDAT_BOOST_SHA1})
    else ()
        set(BoostVersion 1.67.0)
        set(BoostSHA1 694ae3f4f899d1a80eb7a3b31b33be73c423c1ae)
    endif ()
    message(STATUS "Downloading Boost Version ${BoostVersion}")
    include(cmake_modules/add_boost.cmake)
    message(STATUS "Downloading Boost Version ${BoostVersion} - done")
    #    set(BOOST_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/boost)
    #    set(Boost_LIBRARY ${CMAKE_CURRENT_SOURCE_DIR}/boost/stage/lib)
    #    set(Boost_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/boost)
    #    set(Boost_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/boost)
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

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")

#if (NOT TUDAT_DOWNLOAD_AND_BUILD_BOOST)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
#    set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -D_GLIBCXX_USE_CXX11_ABI=0")
#endif ()
include(TudatFindBoost)
message(STATUS "Checking for _GLIBCXX_USE_CXX11_ABI definition...")
if (DEFINED _GLIBCXX_USE_CXX11_ABI)
    if (_GLIBCXX_USE_CXX11_ABI)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=1")
        set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -D_GLIBCXX_USE_CXX11_ABI=1")
        message(STATUS "-D_GLIBCXX_USE_CXX11_ABI=1 added to CMAKE_CXX_FLAGS!")
    elseif (NOT _GLIBCXX_USE_CXX11_ABI)
        message(Boost_VERSION_MAJOR: ${Boost_VERSION_MAJOR})
        message(Boost_VERSION_MINOR: ${Boost_VERSION_MINOR})
        message(TUDAT_BUILD_GNU: ${TUDAT_BUILD_GNU})
        if ((${Boost_VERSION_MAJOR} LESS_EQUAL "1") AND (${Boost_VERSION_MINOR} LESS_EQUAL "67") AND ${TUDAT_BUILD_GNU})
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
            set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -D_GLIBCXX_USE_CXX11_ABI=0")
            message(STATUS "-D_GLIBCXX_USE_CXX11_ABI=0 added to CMAKE_CXX_FLAGS!")
        elseif ((${Boost_VERSION_MAJOR} GREATER_EQUAL "1") AND (${Boost_VERSION_MINOR} GREATER_EQUAL "70") AND ${TUDAT_BUILD_GNU})
            #            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
            #            set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -D_GLIBCXX_USE_CXX11_ABI=0")
            #            message(STATUS "-D_GLIBCXX_USE_CXX11_ABI=0 added to CMAKE_CXX_FLAGS!")
        elseif (${TUDAT_BUILD_MSVC})
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
            set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -D_GLIBCXX_USE_CXX11_ABI=0")
            message(STATUS "-D_GLIBCXX_USE_CXX11_ABI=0 added to CMAKE_CXX_FLAGS!")
        else ()
            message(STATUS "_GLIBCXX_USE_CXX11_ABI was set as 0, but since the boost version is >=1.68|>1.70 [${Boost_VERSION}] and compiler_id is GNU, it's being ignored!")
            add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)
        endif ()
    endif ()

else ()
    message(STATUS "_GLIBCXX_USE_CXX11_ABI was not found.")
endif ()

#if (NOT TUDAT_DOWNLOAD_AND_BUILD_BOOST)
#    if (Boost_MAJOR LESS_EQUAL 1)
#        if (Boost_MINOR LESS 67)
#            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
#            set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -D_GLIBCXX_USE_CXX11_ABI=0")
#        endif ()
#    endif ()
#endif ()

# Get tudat resource manager dependency.
find_package(TudatResources REQUIRED)
message(${TudatResources_INCLUDE_DIRS})
include_directories(TudatResources_INCLUDE_DIRS)

# Find Eigen3 library on local system and include directories.
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")
if (APPLE OR (WIN32 AND CLANG))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
endif ()


# Sofa dependency if in build settings.
if (TUDAT_BUILD_WITH_PAGMO)
    #
    add_definitions(-DTUDAT_BUILD_WITH_PAGMO1)
    if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
        find_package(Pagmo REQUIRED)
    elseif (NOT ${TUDAT_BUNDLE_BUILD})
        find_package(Pagmo REQUIRED)
    else ()
        set(Pagmo_INCLUDE_DIRS ${Pagmo_BUILD_INCLUDE_DIRS})
    endif ()
else ()
    add_definitions(-DTUDAT_BUILD_WITH_PAGMO=0)
endif ()


# Adds Boost_UNIT_TEST_FRAMEWORK_LIBRARY if required.
# Configure Boost libraries.

#if (TUDAT_DOWNLOAD_AND_BUILD_BOOST)
#    include(TudatFindBoost)
#    #    # Configure Boost libraries.
#    #    set(Boost_USE_STATIC_LIBS ON)
#    #    set(Boost_USE_MULTITHREADED ON)
#    #    set(Boost_USE_STATIC_RUNTIME ON)
#    #    # Find Boost libraries on local system.
#    #    find_package(Boost ${BoostVersion} COMPONENTS ${_TUDAT_REQUIRED_BOOST_LIBS} REQUIRED)
#else ()
#    include(TudatFindBoost)
#endif ()
if (UNIX OR (WIN32 AND CLANG))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif ()

# CSpice dependency
find_package(CSpice REQUIRED 1.0.0)

if (TUDAT_BUILD_WITH_SOFA_INTERFACE)
    find_package(Sofa REQUIRED 1.0.0)
    add_definitions(-DTUDAT_BUILD_WITH_SOFA_INTERFACE=1)
endif ()

find_package(nrlmsise00 REQUIRED 0.1)
message(STATUS ${NRLMSISE00_LIBRARIES})
add_definitions(-DTUDAT_BUILD_WITH_NRLMSISE=1)

if (TUDAT_BUILD_WITH_FFTW3)
    if(NOT APPLE)
        find_package(FFTW3)
        message(STATUS "FFTW " ${FFTW3_LIBRARIES})
        include_directories(SYSTEM AFTER "${FFTW3_INCLUDE_DIRS}")
        if(WIN32)
            set (FFTW3_LIBRARIES "${FFTW3_LIBRARY_DIRS}/fftw3.lib")
        else()
            set (FFTW3_LIBRARIES "${FFTW3_LIBRARY_DIRS}/libfftw3.so")
        endif()
        add_definitions(-DTUDAT_BUILD_WITH_FFTW3=1)
    else( )
        find_package(FFTW3)
        message(STATUS ${FFTW3_LIBRARIES})
        include_directories(SYSTEM AFTER "${FFTW3_INCLUDE_DIRS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${FFTW3_INCLUDE_DIR}\"")
        set (FFTW3_LIBRARIES "${FFTW3_LIBRARY_DIRS}/libfftw3.dylib")
        add_definitions(-DTUDAT_BUILD_WITH_FFTW3=1)
    endif( )
else( )
    add_definitions(-DTUDAT_BUILD_WITH_FFTW3=0)
endif( )

# JSON
if (TUDAT_BUILD_WITH_JSON_INTERFACE)
    find_package(nlohmann_json REQUIRED 3.7.3)
endif ()

# PRECOMPILE EXTENDED PRECISION STUFF
if (NOT TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS)
    message(STATUS "Extended precision propagation disabled!")
    add_definitions(-DTUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=0)
else ()
    message(STATUS "Extended precision propagation enabled!")
    add_definitions(-DTUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=1)
endif ()

if (NOT TUDAT_BUILD_WITH_ESTIMATION_TOOLS)
    add_definitions(-DTUDAT_BUILD_WITH_ESTIMATION_TOOLS=0)
else ()
    add_definitions(-DTUDAT_BUILD_WITH_ESTIMATION_TOOLS=1)
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
add_subdirectory(src/utils)
#add_subdirectory(examples)
#add_subdirectory(applications)


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

#    Copyright (c) 2010-2018, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.

# Specify minimum CMake version required to compile Tudat.
cmake_minimum_required( VERSION 2.6 )

# Specific Tudat project name.
project( Tudat )

# Load UserSettings.txt
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
 set(BUILD_STYLE "standalone")
 include("${CMAKE_CURRENT_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
else()
 set(BUILD_STYLE "part of ${CMAKE_PROJECT_NAME}")
 include("${CMAKE_CURRENT_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
 include("${CMAKE_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
 STRING(REGEX REPLACE ${CMAKE_SOURCE_DIR} "" RELATIVE_PROJECT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
 set(RELATIVE_PROJECT_PATH "${RELATIVE_PROJECT_PATH}" CACHE STRING "Relative path wrt to project for function")
 # message(STATUS "Relative path (wrt to project): ${RELATIVE_PROJECT_PATH}")
endif()

# Set CMake build-type. If it not supplied by the user (either directly as an argument of through
# the "UserSettings.txt" file, the default built type is "Release".
if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
 set(CMAKE_BUILD_TYPE Release)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
 set(CMAKE_BUILD_TYPE Debug)
endif()

message(STATUS "<< ${PROJECT_NAME} (${CMAKE_BUILD_TYPE} - ${BUILD_STYLE}) >>")

# Add local module path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/External/CMake")
message(STATUS "CMake Module path(s): ${CMAKE_MODULE_PATH}")

# Set compiler based on preferences (e.g. USE_CLANG) and system.
include(compiler)

# Set root-directory for code to current source directory.
set(CODEROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")

# Define the basic input (trunk).
set(SRCROOT "${CODEROOT}/Tudat")

# Set lib and bin directories where static libraries and unit tests are built.
set(LIBROOT "${CODEROOT}/lib")
set(BINROOT "${CODEROOT}/bin")

# Set testing options based on platform.
enable_testing( )

# Set the global macros for setting up targets.
macro(setup_tudat_library_target target_name CUSTOM_OUTPUT_PATH)
 set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIBROOT}")
 set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIBROOT}")
endmacro(setup_tudat_library_target)

macro(setup_custom_test_program target_name CUSTOM_OUTPUT_PATH)
 set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/unit_tests")
 get_property(CUSTOM_TEST_PROGRAM_NAME TARGET ${target_name} PROPERTY OUTPUT_NAME)
 add_test("${target_name}" "${BINROOT}/unit_tests/${target_name}")
endmacro(setup_custom_test_program)

# Set the main sub-directories.
set(ASTRODYNAMICSDIR "/Astrodynamics")
set(BASICSDIR "/Basics")
set(EXTERNALDIR "/External")
set(INPUTOUTPUTDIR "/InputOutput")
set(MATHEMATICSDIR "/Mathematics")
set(SIMULATIONSETUPDIR "/SimulationSetup")


# Include the top-level directories.
include_directories(AFTER
 "${CODEROOT}"
)

# Find Eigen3 library on local system.
find_package(Eigen3 REQUIRED)

# Include Eigen3 directories.
# Set CMake flag to suppress Eigen warnings (platform-dependent solution).
if(NOT APPLE OR APPLE_INCLUDE_FORCE)
 include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")
else( )
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
endif( )

# Configure Boost libraries.
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)

# Find Boost libraries on local system.
find_package(Boost 1.45.0 COMPONENTS date_time system unit_test_framework filesystem regex REQUIRED)

# Include Boost directories.
# Set CMake flag to suppress Boost warnings (platform-dependent solution).
if(NOT APPLE OR APPLE_INCLUDE_FORCE)
 include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")
else( )
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif( )

# Add an option to toggle the generation of the API documentation.
# If documentation should be built, find Doxygen package and setup config file.
option(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)
if(BUILD_DOCUMENTATION)
 find_package(Doxygen)

 if(NOT DOXYGEN_FOUND)
   message(FATAL_ERROR
           "Doxygen is needed to build the documentation. Please install it correctly.")
 endif()

 # Configure the Template Doxyfile for our specific project.
 configure_file(${CODEROOT}/Doxyfile.in ${PROJECT_BINARY_DIR}/Doxyfile @ONLY IMMEDIATE)

 # Add a custom target to run Doxygen when ever the project is built.
 # IF you do NOT want the documentation to be generated EVERY time you build the project
 # then leave out the 'ALL' keyword from the above command.
 add_custom_target (docs ALL COMMAND ${DOXYGEN_EXECUTABLE}
                    ${PROJECT_BINARY_DIR}/Doxyfile
                    SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)

 # Install documentation.
 install(DIRECTORY "${CODEROOT}/Documentation"
         DESTINATION "${TUDAT_DISTRIBUTION_PATH}")
endif()

#
# CSPICE
#
# Set whether to use the SPICE library integration with Tudat or not. If it not supplied by the
# user (either directly as an argument or through the "UserSettings.txt" file, the default setting
# is "OFF". NB: It is important this section is placed AFTER the "UserSettings.txt" file is
# included, to ensure that that the correct indiciation is given for whether SPICE should be
# included or not has been selected.
option(USE_CSPICE "build Tudat with CSPICE enabled" ON) # OFF is the default
if(NOT USE_CSPICE)
 message(STATUS "SPICE disabled!")
 add_definitions(-DUSE_CSPICE=0)
else()
 message(STATUS "SPICE enabled!")
add_definitions(-DUSE_CSPICE=1)

 # Find CSPICE library on local system.
 find_package(Spice)

 # Include CSpice directories.
 if(NOT APPLE OR APPLE_INCLUDE_FORCE)
   include_directories(SYSTEM AFTER "${SPICE_INCLUDE_DIR}")
 else( )
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${SPICE_INCLUDE_DIR}\"")
 endif( )
endif()


#
# JSON
#
# Set whether to use the JSON library.
option(USE_JSON "build Tudat with JSON enabled" ON)
if(NOT USE_JSON)
 message(STATUS "JSON disabled!")
 add_definitions(-DUSE_JSON=0)
else()
 message(STATUS "JSON enabled!")
 add_definitions(-DUSE_JSON=1)

 # Find JSON library on local system.
 find_package(nlohmann_json)

 # Include JSON directories.
 if(NOT APPLE OR APPLE_INCLUDE_FORCE)
   include_directories(SYSTEM AFTER "${JSON_INCLUDE_DIR}")
 else( )
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${JSON_INCLUDE_DIR}\"")
 endif( )
endif()


#
# SOFA
#
# Set whether to use the SOFA library integration with Tudat or not. If it not supplied by the
# user (either directly as an argument or through the "UserSettings.txt" file, the default setting
# is "OFF". NB: It is important this section is placed AFTER the "UserSettings.txt" file is
# included, to ensure that that the correct indiciation is given for whether SOFA should be
# included or not has been selected.
option(USE_SOFA "build Tudat with Sofa enabled" ON) # OFF is the default
if(NOT USE_SOFA)
 message(STATUS "SOFA disabled!")
 add_definitions(-DUSE_SOFA=0)
else()
 message(STATUS "SOFA enabled!")
add_definitions(-DUSE_SOFA=1)

 # Find SOFA library on local system.
 find_package(Sofa)

 # Include Sofa directories.
 if(NOT APPLE OR APPLE_INCLUDE_FORCE)
   include_directories(SYSTEM AFTER "${SOFA_INCLUDE_DIR}")
 else( )
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${SOFA_INCLUDE_DIR}\"")
 endif( )
endif()

#
# NRLMSISE-00
#
# Set whether to use the NRLMSISE-00 library integration with Tudat or not. If it not supplied by the
# user (either directly as an argument or through the "UserSettings.txt" file, the default setting
# is "ON").
option(USE_NRLMSISE00 "build Tudat with NRLMSISE-00 enabled" ON)
if(NOT USE_NRLMSISE00)
 message(STATUS "NRLMSISE-00 disabled!")
 add_definitions(-DUSE_NRLMSISE00=0)
else()
 message(STATUS "NRLMSISE-00 enabled!")
 add_definitions(-DUSE_NRLMSISE00=1)
 # Find USE_NRLMSISE00 library on local system.
 find_package(NRLMSISE00)

 # Include NRLMSISE00 directories.
 if(NOT APPLE OR APPLE_INCLUDE_FORCE)
   include_directories(SYSTEM AFTER "${NRLMSISE00_INCLUDE_DIR}")
 else( )
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${NRLMSISE00_INCLUDE_DIR}\"")
 endif( )
endif()

#
# GSL
#
# Set whether to use the GSL library integration with Tudat or not. If it not supplied by the
# user (either directly as an argument or through the "UserSettings.txt" file, the default setting
# is "OFF").
option(USE_GSL "build Tudat with GSL enabled" OFF)
if(NOT USE_GSL)
 message(STATUS "GSL disabled!")
 add_definitions(-DUSE_GSL=0)
else()
 message(STATUS "GSL enabled!")
 add_definitions(-DUSE_GSL=1)
 # Find USE_NRLMSISE00 library on local system.
 find_package(GSL)

 # Include NRLMSISE00 directories.
 if(NOT APPLE OR APPLE_INCLUDE_FORCE)
   include_directories(SYSTEM AFTER "${GSL_INCLUDE_DIR}")
 else( )
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${GSL_INCLUDE_DIR}\"")
 endif( )
endif()

option(COMPILE_HIGH_ACCURACY_ESTIMATION_TESTS  "Compiling unit tests for state estimation. These may cause excessive (>3 GB)) RAM usage with gcc/mingw." ON)
option(COMPILE_PROPAGATION_TESTS "Compiling unit tests involving long (> 30 s) propagations. Total unit test run time may be > 5-10 minutes." ON)

# Create lists of static libraries for ease of use
list(APPEND TUDAT_EXTERNAL_LIBRARIES "")
list(APPEND TUDAT_EXTERNAL_INTERFACE_LIBRARIES "")
list(APPEND TUDAT_ITRS_LIBRARIES "")

if(USE_SOFA)
 list(APPEND TUDAT_EXTERNAL_LIBRARIES sofa)
 list(APPEND TUDAT_EXTERNAL_INTERFACE_LIBRARIES tudat_sofa_interface )
 list(APPEND TUDAT_ITRS_LIBRARIES tudat_earth_orientation )
endif()

if(USE_CSPICE)
 list(APPEND TUDAT_EXTERNAL_LIBRARIES cspice)
 list(APPEND TUDAT_EXTERNAL_INTERFACE_LIBRARIES tudat_spice_interface )
endif()

if(USE_JSON)
  list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES nlohmann_json)
endif()

if(USE_NRLMSISE00)
 list(APPEND TUDAT_EXTERNAL_LIBRARIES nrlmsise00)
endif()

if(USE_GSL)
 list(APPEND TUDAT_EXTERNAL_LIBRARIES gsl)
endif()


list(APPEND TUDAT_PROPAGATION_LIBRARIES tudat_simulation_setup tudat_ground_stations tudat_propagators
    tudat_aerodynamics tudat_system_models tudat_geometric_shapes tudat_relativity tudat_gravitation tudat_mission_segments
    tudat_electro_magnetism tudat_propulsion tudat_ephemerides ${TUDAT_ITRS_LIBRARIES} tudat_numerical_integrators tudat_reference_frames
     tudat_statistics tudat_propagators ${TUDAT_EXTERNAL_INTERFACE_LIBRARIES} tudat_basic_astrodynamics tudat_root_finders tudat_basic_mathematics
     tudat_input_output tudat_basics ${TUDAT_EXTERNAL_LIBRARIES})

list(APPEND TUDAT_ESTIMATION_LIBRARIES tudat_simulation_setup tudat_observation_models tudat_ground_stations tudat_acceleration_partials
    tudat_observation_partials tudat_estimatable_parameters tudat_orbit_determination  tudat_propagators
    tudat_aerodynamics tudat_system_models tudat_geometric_shapes tudat_relativity tudat_gravitation tudat_mission_segments
    tudat_electro_magnetism tudat_propulsion tudat_ephemerides ${TUDAT_ITRS_LIBRARIES} tudat_numerical_integrators tudat_reference_frames
    tudat_statistics tudat_propagators ${TUDAT_EXTERNAL_INTERFACE_LIBRARIES} tudat_basic_astrodynamics tudat_root_finders tudat_basic_mathematics
    tudat_input_output tudat_basics ${TUDAT_EXTERNAL_LIBRARIES})



# Set sub-directories.
set(SUBDIRS ${APPLICATIONSDIR} ${ASTRODYNAMICSDIR} ${BASICSDIR} ${INPUTOUTPUTDIR} ${MATHEMATICSDIR} ${EXTERNALDIR} ${SIMULATIONSETUPDIR})

if(USE_JSON)
  # Set JSON interface directory.
  set(JSONINTERFACEDIR "/JsonInterface")

  # Add subdirectories.
  list(APPEND SUBDIRS ${JSONINTERFACEDIR})
endif()

# Add sub-directories to CMake process.
foreach(CURRENT_SUBDIR ${SUBDIRS})
  add_subdirectory("${SRCROOT}${CURRENT_SUBDIR}")
endforeach(CURRENT_SUBDIR)

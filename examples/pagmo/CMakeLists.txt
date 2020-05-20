 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

# Specify minimum cmake version required.
cmake_minimum_required(VERSION 2.6)

# Specify project name.
project(SatellitePropagatorExamples)

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

# Set cmake build-type. If it not supplied by the user (either directly as an argument of through
# the "UserSettings.txt" file, the default built type is "Release".
if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
  set(CMAKE_BUILD_TYPE Release)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_BUILD_TYPE Debug)
endif()

message(STATUS "<< ${PROJECT_NAME} (${CMAKE_BUILD_TYPE} - ${BUILD_STYLE}) >>")

# Add local module path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMakeModules")
message(STATUS "CMake Module path(s): ${CMAKE_MODULE_PATH}")

# Set compiler based on preferences (e.g. USE_CLANG) and system.
include(compiler)

# Define the directory with the source code.
set(SRCROOT "${CMAKE_CURRENT_SOURCE_DIR}")

# Define the code root directory.
set(CODEROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")

# Set testing options based on platform.
enable_testing()

# Set lib and bin directories where static libraries and unit tests are built.
if(NOT LIBROOT)
  set(LIBROOT "${CODEROOT}/lib")
endif()
if(NOT BINROOT)
  set(BINROOT "${CODEROOT}/bin")
endif()

# Set the global macros for setting up targets.
macro(setup_executable_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/applications")
  install(TARGETS ${target_name} RUNTIME DESTINATION "${BINROOT}/applications")
endmacro(setup_executable_target)

macro(setup_library_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIBROOT}")
  set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIBROOT}")
endmacro(setup_library_target)

macro(setup_unit_test_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/unit_tests")
  get_property(CUSTOM_TEST_PROGRAM_NAME TARGET ${target_name} PROPERTY OUTPUT_NAME)
  add_test("${target_name}" "${BINROOT}/unit_tests/${target_name}")
endmacro(setup_unit_test_target)

# Include the top-level directories.
include_directories(AFTER
  "${CODEROOT}"
)

# Find Eigen3 library on local system.
find_package(Eigen3 REQUIRED)

# Include Eigen3 directories.
# Set cmake flag to suppress Eigen warnings (platform-dependent solution).
if(NOT APPLE OR APPLE_INCLUDE_FORCE)
  include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
endif()

# Configure Boost libraries.
if(NOT Boost_USE_STATIC_LIBS)
  set(Boost_USE_STATIC_LIBS ON)
endif()
if(NOT Boost_USE_MULTITHREADED)
  set(Boost_USE_MULTITHREADED ON)
endif()
if(NOT Boost_USE_STATIC_RUNTIME)
  set(Boost_USE_STATIC_RUNTIME ON)
endif()

# Find Boost libraries on local system.
find_package(Boost 1.55.0
             COMPONENTS thread date_time system unit_test_framework filesystem regex REQUIRED)

# Include Boost directories.
# Set cmake flag to suppress Boost warnings (platform-dependent solution).
if(NOT APPLE OR APPLE_INCLUDE_FORCE)
  include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif()

# Find Tudat library on local system.
find_package(Tudat 2.0 REQUIRED)

# Include Tudat directories.
# Set cmake flag to suppress Tudat warnings (platform-dependent solution).
if(NOT APPLE OR APPLE_INCLUDE_FORCE)
  include_directories(SYSTEM AFTER "${TUDAT_INCLUDE_DIR}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${TUDAT_INCLUDE_DIR}\"")
endif()

  # Find CSPICE library on local system.
  find_package(Spice)

  # Include CSpice directories.
  if(NOT APPLE OR APPLE_INCLUDE_FORCE)
    include_directories(SYSTEM AFTER "${SPICE_INCLUDE_DIR}")
  else( )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${SPICE_INCLUDE_DIR}\"")
  endif( )

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
endif( )

# Find PaGMO library on local system.
if( USE_PAGMO )
  find_package(PaGMO2)

  # Include PaGMO directories.
  if(NOT APPLE OR APPLE_INCLUDE_FORCE)
    include_directories(SYSTEM AFTER "${PAGMO_INCLUDE_DIR}")
  else( )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${PAGMO_INCLUDE_DIR}\"")
  endif( )
endif( )

# Set compiler based on preferences (e.g. USE_CLANG) and system.
#include(tudatLinkLibraries)

## Add Tudat pagmo examples
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/Problems")

# Disable the CEC2013 example application on Windows as CEC2013 is not supported by MinGW
if( NOT MINGW )
    add_executable(application_PagmoCEC2013Comparison "${SRCROOT}/cec2013OptimizerComparison.cpp")
    setup_executable_target(application_PagmoCEC2013Comparison "${SRCROOT}")
    target_link_libraries(application_PagmoCEC2013Comparison pagmo pagmo2_library_example_problems tudat_input_output pthread ${Boost_LIBRARIES} )
endif()

add_executable(application_PagmoHimmelblauOptimization "${SRCROOT}/himmelblauOptimization.cpp")
setup_executable_target(application_PagmoHimmelblauOptimization "${SRCROOT}")
target_link_libraries(application_PagmoHimmelblauOptimization pagmo pagmo2_library_example_problems tudat_input_output pthread ${Boost_LIBRARIES} )

add_executable(application_PagmoEarthMarsTransferExample "${SRCROOT}/earthMarsTransferExample.cpp")
setup_executable_target(application_PagmoEarthMarsTransferExample "${SRCROOT}")
target_link_libraries(application_PagmoEarthMarsTransferExample pagmo pagmo2_library_example_problems ${Tudat_PROPAGATION_LIBRARIES}
     pthread ${Boost_LIBRARIES} )

 add_executable(application_PagmoMgaTransferExample "${SRCROOT}/mgaTransferExample.cpp")
 setup_executable_target(application_PagmoMgaTransferExample "${SRCROOT}")
 target_link_libraries(application_PagmoMgaTransferExample pagmo pagmo2_library_example_problems ${Tudat_PROPAGATION_LIBRARIES}
      pthread ${Boost_LIBRARIES} )


 add_executable(application_PagmoMultiObjectiveEarthMarsTransferExample "${SRCROOT}/multiObjectiveEarthMarsTransferExample.cpp")
 setup_executable_target(application_PagmoMultiObjectiveEarthMarsTransferExample "${SRCROOT}")
 target_link_libraries(application_PagmoMultiObjectiveEarthMarsTransferExample pagmo pagmo2_library_example_problems ${Tudat_PROPAGATION_LIBRARIES}
      pthread ${Boost_LIBRARIES} )

add_executable(application_PagmoPropagationTargetingExample "${SRCROOT}/propagationTargetingExample.cpp")
setup_executable_target(application_PagmoPropagationTargetingExample "${SRCROOT}")
target_link_libraries(application_PagmoPropagationTargetingExample pagmo pagmo2_library_example_problems ${Tudat_PROPAGATION_LIBRARIES}
      pthread ${Boost_LIBRARIES} )

if( BUILD_WITH_SOCIS_2019 )

 add_executable(application_PagmoHodographicShapingTrajectoryExample "${SRCROOT}/hodographicShapingTrajectoryExample.cpp")
 setup_executable_target(application_PagmoHodographicShapingTrajectoryExample "${SRCROOT}")
 target_link_libraries(application_PagmoHodographicShapingTrajectoryExample tudat_shape_based_methods tudat_low_thrust_trajectories
     tudat_numerical_quadrature pagmo ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

 add_executable(application_PagmoHodographicShapingFullOptimisation "${SRCROOT}/hodographicShapingFullOptimisationExample.cpp")
 setup_executable_target(application_PagmoHodographicShapingFullOptimisation "${SRCROOT}")
 target_link_libraries(application_PagmoHodographicShapingFullOptimisation tudat_shape_based_methods tudat_low_thrust_trajectories
     tudat_numerical_quadrature pagmo ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

# # Add low-thrust trajectory design example
# add_executable(application_PagmoSimsFlanaganTrajectoryExample "${SRCROOT}/simsFlanaganTrajectoryExample.cpp")
# setup_executable_target(application_PagmoSimsFlanaganTrajectoryExample "${SRCROOT}")
# target_link_libraries(application_PagmoSimsFlanaganTrajectoryExample pagmo2_library_example_problems tudat_low_thrust_trajectories tudat_shape_based_methods
#     tudat_numerical_quadrature ${Tudat_PROPAGATION_LIBRARIES} pthread pagmo ${Boost_LIBRARIES} )

endif( )

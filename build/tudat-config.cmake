# defined since 2.8.3
if (CMAKE_VERSION VERSION_LESS 2.8.3)
  get_filename_component (CMAKE_CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_FILE} PATH)
endif ()

# Temporarily modify CMAKE_MODULE_PATH for cmake files in current dir.
set(_TUDAT_CONFIG_OLD_MODULE_PATH "${CMAKE_MODULE_PATH}")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

# Find dependencies.
include(CMakeFindDependencyMacro)
find_dependency(CSpice)
find_dependency(Sofa)
#find_dependency(Eigen3)
#efind_dependency(Boost)
#set(_TUDAT_FIND_BOOST_UNIT_TEST_FRAMEWORK ON)
#include(TudatFindBoost)



# Tell the user project where to find our headers and libraries
set (Tudat_VERSION "")
set (Tudat_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/../../include")
set (Tudat_LIBRARY_DIRS "${CMAKE_CURRENT_LIST_DIR}/../")
set (Tudat_DATA_DIRS "${CMAKE_CURRENT_LIST_DIR}/../../data")

# Configure file path for tudat data loading.
##configure_file(
##        "${CMAKE_CURRENT_LIST_DIR}/paths.hpp.in"
##        "${Tudat_INCLUDE_DIRS}/tudat/paths.hpp" @ONLY
##)

# List of compilation flags -DTOTO to export
set (Tudat_DEFINITIONS " -DTUDAT_BUILD_WITH_FILTERS=OFF -DTUDAT_BUILD_WITH_SOFA_INTERFACE=ON -DTUDAT_BUILD_WITH_FFTW3=OFF -DTUDAT_BUILD_WITH_JSON_INTERFACE=OFF -DTUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=ON")

# Optional dependencies.


# Allows loading CSpice settings from another project
set (Tudat_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")

# Our library dependencies (contains definitions for IMPORTED targets)
include ("${CMAKE_CURRENT_LIST_DIR}/tudat_export.cmake")

# These are IMPORTED targets created by tudat_targets.cmake
set (Tudat_PROPAGATION_LIBRARIES "Tudat::tudat_propagation_setup;Tudat::tudat_shape_based_methods;Tudat::tudat_low_thrust_trajectories;Tudat::tudat_environment_setup;Tudat::tudat_ground_stations;Tudat::tudat_aerodynamics;Tudat::tudat_system_models;Tudat::tudat_geometric_shapes;Tudat::tudat_relativity;Tudat::tudat_gravitation;Tudat::tudat_mission_segments;Tudat::tudat_electromagnetism;Tudat::tudat_propulsion;Tudat::tudat_ephemerides;Tudat::tudat_earth_orientation;Tudat::tudat_numerical_integrators;Tudat::tudat_reference_frames;Tudat::tudat_statistics;Tudat::tudat_propagators;Tudat::tudat_spice_interface;Tudat::tudat_sofa_interface;NRLMSISE00::nrlmsise00;Tudat::tudat_basic_astrodynamics;Tudat::tudat_numerical_quadrature;Tudat::tudat_interpolators;Tudat::tudat_root_finders;Tudat::tudat_basic_mathematics;Tudat::tudat_input_output;Tudat::tudat_basics;Tudat::tudat_data")
set (Tudat_ESTIMATION_LIBRARIES "Tudat::tudat_estimation_setup;Tudat::tudat_propagation_setup;Tudat::tudat_shape_based_methods;Tudat::tudat_low_thrust_trajectories;Tudat::tudat_environment_setup;Tudat::tudat_observation_models;Tudat::tudat_ground_stations;Tudat::tudat_acceleration_partials;Tudat::tudat_torque_partials;Tudat::tudat_observation_partials;Tudat::tudat_orbit_determination;Tudat::tudat_estimatable_parameters;Tudat::tudat_aerodynamics;Tudat::tudat_system_models;Tudat::tudat_geometric_shapes;Tudat::tudat_relativity;Tudat::tudat_gravitation;Tudat::tudat_mission_segments;Tudat::tudat_electromagnetism;Tudat::tudat_propulsion;Tudat::tudat_ephemerides;Tudat::tudat_earth_orientation;Tudat::tudat_numerical_integrators;Tudat::tudat_reference_frames;Tudat::tudat_statistics;Tudat::tudat_propagators;Tudat::tudat_spice_interface;Tudat::tudat_sofa_interface;NRLMSISE00::nrlmsise00;Tudat::tudat_basic_astrodynamics;Tudat::tudat_numerical_quadrature;Tudat::tudat_interpolators;Tudat::tudat_root_finders;Tudat::tudat_basic_mathematics;Tudat::tudat_input_output;Tudat::tudat_basics;Tudat::tudat_data")

if (CMAKE_VERSION VERSION_LESS 2.8.3)
  set (CMAKE_CURRENT_LIST_DIR)
endif ()

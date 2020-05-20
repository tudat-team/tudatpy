
# defined since 2.8.3
if (CMAKE_VERSION VERSION_LESS 2.8.3)
  get_filename_component (CMAKE_CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_FILE} PATH)
endif ()

# Get dependencies using CMakeFindDependencyMacro.

set(_TUDAT_FIND_BOOST_UNIT_TEST_FRAMEWORK ON)
include(CMakeFindDependencyMacro)
find_dependency(CSpice)
find_dependency(Sofa)
find_dependency(Eigen3)
include(TudatFindBoost)

# Tell the user project where to find our headers and libraries
set (Tudat_VERSION "4.0.0")
set (Tudat_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/")
set (Tudat_LIBRARY_DIRS "${CMAKE_CURRENT_LIST_DIR}/")
set (Tudat_DATA_DIRS "${CMAKE_CURRENT_LIST_DIR}/")

# List of compilation flags -DTOTO to export
set (Tudat_DEFINITIONS "")

# Optional dependencies.


# Allows loading CSpice settings from another project
# set (Tudat_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")

# Our library dependencies (contains definitions for IMPORTED targets)
include ("${CMAKE_CURRENT_LIST_DIR}/tudat_targets.cmake")

# These are IMPORTED targets created by tudat_targets.cmake
set (Tudat_PROPAGATION_LIBRARIES "tudat_propagation_setup;tudat_trajectory_design;tudat_environment_setup;tudat_ground_stations;tudat_propagators;tudat_aerodynamics;tudat_system_models;tudat_geometric_shapes;tudat_relativity;tudat_gravitation;tudat_mission_segments;tudat_electromagnetism;tudat_propulsion;tudat_ephemerides;tudat_earth_orientation;tudat_numerical_integrators;tudat_reference_frames;tudat_statistics;tudat_propagators;tudat_sofa_interface;tudat_spice_interface;tudat_basic_astrodynamics;tudat_interpolators;tudat_root_finders;tudat_basic_mathematics;tudat_input_output;tudat_basics")
set (Tudat_ESTIMATION_LIBRARIES "tudat_propagation_setup;tudat_trajectory_design;tudat_environment_setup;tudat_ground_stations;tudat_propagators;tudat_aerodynamics;tudat_system_models;tudat_geometric_shapes;tudat_relativity;tudat_gravitation;tudat_mission_segments;tudat_electromagnetism;tudat_propulsion;tudat_ephemerides;tudat_earth_orientation;tudat_numerical_integrators;tudat_reference_frames;tudat_statistics;tudat_propagators;tudat_sofa_interface;tudat_spice_interface;tudat_basic_astrodynamics;tudat_interpolators;tudat_root_finders;tudat_basic_mathematics;tudat_input_output;tudat_basics")

if (CMAKE_VERSION VERSION_LESS 2.8.3)
  set (CMAKE_CURRENT_LIST_DIR)
endif ()


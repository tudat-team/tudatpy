 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #


 set(SHAPEBASEDMETHODSDIR "${LOWTHRUSTTRAJECTORIESDIR}/ShapeBasedMethods")

 # Add subdirectories.
 add_subdirectory("${SRCROOT}${SHAPEBASEDMETHODSDIR}")

 # Get target properties for static libraries.
 get_target_property(SHAPEBASEDMETHOSDSOURCES tudat_shape_based_methods SOURCES)

# Set the source files.
set(LOWTHRUSTTRAJECTORIES_SOURCES
  "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/lowThrustLegSettings.cpp"
  "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/lowThrustLeg.cpp"
  "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/simsFlanaganLeg.cpp"
  "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/simsFlanagan.cpp"
  "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/simsFlanaganOptimisationSetup.cpp"
  "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/hybridMethodLeg.cpp"
  "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/hybridOptimisationSetup.cpp"
  "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/hybridMethod.cpp"
  "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/lowThrustOptimisationSetup.cpp"
#  "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/optimisationSettings.cpp"
)

# Set the header files.
set(LOWTHRUSTTRAJECTORIES_HEADERS
  "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/lowThrustLegSettings.h"
  "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/lowThrustLeg.h"
  "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/simsFlanaganLeg.h"
  "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/simsFlanagan.h"
  "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/simsFlanaganOptimisationSetup.h"
  "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/hybridMethodLeg.h"
  "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/hybridOptimisationSetup.h"
  "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/hybridMethod.h"
  "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/lowThrustOptimisationSetup.h"
#  "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/optimisationSettings.h"
)

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

# Add static libraries.
add_library(tudat_low_thrust_trajectories STATIC ${LOWTHRUSTTRAJECTORIES_SOURCES} ${LOWTHRUSTTRAJECTORIES_HEADERS})
setup_tudat_library_target(tudat_low_thrust_trajectories "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}")

include(tudatLinkLibraries)

## Add unit tests.
add_executable(test_SimsFlanagan "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/UnitTests/unitTestSimsFlanagan.cpp")
setup_custom_test_program(test_SimsFlanagan "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}")
target_link_libraries(test_SimsFlanagan tudat_shape_based_methods tudat_low_thrust_trajectories tudat_numerical_quadrature ${TUDAT_PROPAGATION_LIBRARIES}  ${Boost_LIBRARIES})

add_executable(test_HybridMethod "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}/UnitTests/unitTestHybridMethod.cpp")
setup_custom_test_program(test_HybridMethod "${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}")
target_link_libraries(test_HybridMethod tudat_low_thrust_trajectories tudat_numerical_quadrature ${TUDAT_PROPAGATION_LIBRARIES}  ${Boost_LIBRARIES})

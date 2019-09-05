 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Set the source files.
set(LOWTHRUSTDIRECTMETHODS_SOURCES
  "${SRCROOT}${LOWTHRUSTDIRECTMETHODSDIR}/simsFlanaganLeg.cpp"
  "${SRCROOT}${LOWTHRUSTDIRECTMETHODSDIR}/simsFlanagan.cpp"
  "${SRCROOT}${LOWTHRUSTDIRECTMETHODSDIR}/simsFlanaganOptimisationSetup.cpp"
  "${SRCROOT}${LOWTHRUSTDIRECTMETHODSDIR}/hybridMethodLeg.cpp"
  "${SRCROOT}${LOWTHRUSTDIRECTMETHODSDIR}/hybridOptimisationSetup.cpp"
  "${SRCROOT}${LOWTHRUSTDIRECTMETHODSDIR}/hybridMethod.cpp"
)

# Set the header files.
set(LOWTHRUSTDIRECTMETHODS_HEADERS
  "${SRCROOT}${LOWTHRUSTDIRECTMETHODSDIR}/simsFlanaganLeg.h"
  "${SRCROOT}${LOWTHRUSTDIRECTMETHODSDIR}/simsFlanagan.h"
  "${SRCROOT}${LOWTHRUSTDIRECTMETHODSDIR}/simsFlanaganOptimisationSetup.h"
  "${SRCROOT}${LOWTHRUSTDIRECTMETHODSDIR}/hybridMethodLeg.h"
  "${SRCROOT}${LOWTHRUSTDIRECTMETHODSDIR}/hybridOptimisationSetup.h"
  "${SRCROOT}${LOWTHRUSTDIRECTMETHODSDIR}/hybridMethod.h"
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
add_library(tudat_low_thrust_direct_methods STATIC ${LOWTHRUSTDIRECTMETHODS_SOURCES} ${LOWTHRUSTDIRECTMETHODS_HEADERS})
setup_tudat_library_target(tudat_low_thrust_direct_methods "${SRCROOT}${LOWTHRUSTDIRECTMETHODSDIR}")

include(tudatLinkLibraries)

## Add unit tests.
add_executable(test_SimsFlanagan "${SRCROOT}${LOWTHRUSTDIRECTMETHODSDIR}/UnitTests/unitTestSimsFlanagan.cpp")
setup_custom_test_program(test_SimsFlanagan "${SRCROOT}${LOWTHRUSTDIRECTMETHODSDIR}")
target_link_libraries(test_SimsFlanagan tudat_low_thrust_direct_methods tudat_shape_based_methods tudat_numerical_quadrature ${TUDAT_PROPAGATION_LIBRARIES}  ${Boost_LIBRARIES})

add_executable(test_HybridMethod "${SRCROOT}${LOWTHRUSTDIRECTMETHODSDIR}/UnitTests/unitTestHybridMethod.cpp")
setup_custom_test_program(test_HybridMethod "${SRCROOT}${LOWTHRUSTDIRECTMETHODSDIR}")
target_link_libraries(test_HybridMethod tudat_low_thrust_direct_methods tudat_numerical_quadrature ${TUDAT_PROPAGATION_LIBRARIES}  ${Boost_LIBRARIES})

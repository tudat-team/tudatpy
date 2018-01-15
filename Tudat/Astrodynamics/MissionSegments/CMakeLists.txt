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
set(MISSIONSEGMENTS_SOURCES
  "${SRCROOT}${MISSIONSEGMENTSDIR}/escapeAndCapture.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/gravityAssist.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/improvedInversePolynomialWall.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/lambertTargeterIzzo.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/lambertTargeterGooding.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/lambertRoutines.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/multiRevolutionLambertTargeterIzzo.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/oscillatingFunctionNovak.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/zeroRevolutionLambertTargeterIzzo.cpp"
)

# Set the header files.
set(MISSIONSEGMENTS_HEADERS 
  "${SRCROOT}${MISSIONSEGMENTSDIR}/escapeAndCapture.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/gravityAssist.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/improvedInversePolynomialWall.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/lambertTargeter.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/lambertTargeterIzzo.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/lambertTargeterGooding.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/lambertRoutines.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/multiRevolutionLambertTargeterIzzo.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/oscillatingFunctionNovak.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/zeroRevolutionLambertTargeterIzzo.h"
)

# Add static libraries.
add_library(tudat_mission_segments STATIC ${MISSIONSEGMENTS_SOURCES} ${MISSIONSEGMENTS_HEADERS})
setup_tudat_library_target(tudat_mission_segments "${SRCROOT}${MISSIONSEGMENTSDIR}")

# Add unit tests.
add_executable(test_EscapeAndCapture "${SRCROOT}${MISSIONSEGMENTSDIR}/UnitTests/unitTestEscapeAndCapture.cpp")
setup_custom_test_program(test_EscapeAndCapture "${SRCROOT}${MISSIONSEGMENTSDIR}")
target_link_libraries(test_EscapeAndCapture tudat_mission_segments tudat_gravitation tudat_basic_astrodynamics tudat_basic_mathematics tudat_root_finders tudat_input_output ${Boost_LIBRARIES})

add_executable(test_GravityAssist "${SRCROOT}${MISSIONSEGMENTSDIR}/UnitTests/unitTestGravityAssist.cpp")
setup_custom_test_program(test_GravityAssist "${SRCROOT}${MISSIONSEGMENTSDIR}")
target_link_libraries(test_GravityAssist tudat_gravitation tudat_mission_segments tudat_basic_astrodynamics tudat_basic_mathematics tudat_root_finders tudat_input_output ${Boost_LIBRARIES})

add_executable(test_LambertTargeterIzzo "${SRCROOT}${MISSIONSEGMENTSDIR}/UnitTests/unitTestLambertTargeterIzzo.cpp")
setup_custom_test_program(test_LambertTargeterIzzo "${SRCROOT}${MISSIONSEGMENTSDIR}")
target_link_libraries(test_LambertTargeterIzzo tudat_mission_segments tudat_root_finders tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_LambertTargeterGooding "${SRCROOT}${MISSIONSEGMENTSDIR}/UnitTests/unitTestLambertTargeterGooding.cpp")
setup_custom_test_program(test_LambertTargeterGooding "${SRCROOT}${MISSIONSEGMENTSDIR}")
target_link_libraries(test_LambertTargeterGooding tudat_mission_segments tudat_root_finders tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_LambertTargeter "${SRCROOT}${MISSIONSEGMENTSDIR}/UnitTests/unitTestLambertTargeter.cpp")
setup_custom_test_program(test_LambertTargeter "${SRCROOT}${MISSIONSEGMENTSDIR}")
target_link_libraries(test_LambertTargeter tudat_mission_segments tudat_root_finders tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_LambertRoutines "${SRCROOT}${MISSIONSEGMENTSDIR}/UnitTests/unitTestLambertRoutines.cpp")
setup_custom_test_program(test_LambertRoutines "${SRCROOT}${MISSIONSEGMENTSDIR}")
target_link_libraries(test_LambertRoutines tudat_mission_segments tudat_root_finders tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_ZeroRevolutionLambertTargeterIzzo "${SRCROOT}${MISSIONSEGMENTSDIR}/UnitTests/unitTestZeroRevolutionLambertTargeterIzzo.cpp")
setup_custom_test_program(test_ZeroRevolutionLambertTargeterIzzo "${SRCROOT}${MISSIONSEGMENTSDIR}")
target_link_libraries(test_ZeroRevolutionLambertTargeterIzzo tudat_mission_segments tudat_root_finders tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_MultiRevolutionLambertTargeterIzzo "${SRCROOT}${MISSIONSEGMENTSDIR}/UnitTests/unitTestMultiRevolutionLambertTargeterIzzo.cpp")
setup_custom_test_program(test_MultiRevolutionLambertTargeterIzzo "${SRCROOT}${MISSIONSEGMENTSDIR}")
target_link_libraries(test_MultiRevolutionLambertTargeterIzzo tudat_mission_segments tudat_root_finders tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_MathematicalShapeFunctions "${SRCROOT}${MISSIONSEGMENTSDIR}/UnitTests/unitTestMathematicalShapeFunctions.cpp")
setup_custom_test_program(test_MathematicalShapeFunctions "${SRCROOT}${MISSIONSEGMENTSDIR}")
target_link_libraries(test_MathematicalShapeFunctions tudat_mission_segments tudat_basic_mathematics ${Boost_LIBRARIES})

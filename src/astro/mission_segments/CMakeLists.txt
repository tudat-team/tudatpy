#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Set the source files.
set(mission_segments_SOURCES
        "escapeAndCapture.cpp"
        "gravityAssist.cpp"
        "improvedInversePolynomialWall.cpp"
        "lambertTargeterIzzo.cpp"
        "lambertTargeterGooding.cpp"
        "lambertRoutines.cpp"
        "multiRevolutionLambertTargeterIzzo.cpp"
        "oscillatingFunctionNovak.cpp"
        "zeroRevolutionLambertTargeterIzzo.cpp"
        )

# Set the header files.
set(mission_segments_HEADERS
        "escapeAndCapture.h"
        "gravityAssist.h"
        "improvedInversePolynomialWall.h"
        "lambertTargeter.h"
        "lambertTargeterIzzo.h"
        "lambertTargeterGooding.h"
        "lambertRoutines.h"
        "multiRevolutionLambertTargeterIzzo.h"
        "oscillatingFunctionNovak.h"
        "zeroRevolutionLambertTargeterIzzo.h"
        )

# Add library.
TUDAT_ADD_LIBRARY("mission_segments"
        "${mission_segments_SOURCES}"
        "${mission_segments_HEADERS}"
        PRIVATE_LINKS "${Boost_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )

# Add static libraries.
#add_library(tudat_mission_segments STATIC ${MISSIONSEGMENTS_SOURCES} ${MISSIONSEGMENTS_HEADERS})
#setup_tudat_library_target(tudat_mission_segments "${SRCROOT}${MISSIONSEGMENTSDIR}")

# Add unit tests.
#add_executable(test_EscapeAndCapture "${SRCROOT}${MISSIONSEGMENTSDIR}/tests/unitTestEscapeAndCapture.cpp")
#setup_custom_test_program(test_EscapeAndCapture "${SRCROOT}${MISSIONSEGMENTSDIR}")
#target_link_libraries(test_EscapeAndCapture tudat_mission_segments tudat_gravitation tudat_basic_astrodynamics tudat_basic_mathematics tudat_root_finders tudat_input_output ${Boost_LIBRARIES})
#
#add_executable(test_GravityAssist "${SRCROOT}${MISSIONSEGMENTSDIR}/tests/unitTestGravityAssist.cpp")
#setup_custom_test_program(test_GravityAssist "${SRCROOT}${MISSIONSEGMENTSDIR}")
#target_link_libraries(test_GravityAssist tudat_gravitation tudat_mission_segments tudat_basic_astrodynamics tudat_basic_mathematics tudat_root_finders tudat_input_output ${Boost_LIBRARIES})
#
#add_executable(test_LambertTargeterIzzo "${SRCROOT}${MISSIONSEGMENTSDIR}/tests/unitTestLambertTargeterIzzo.cpp")
#setup_custom_test_program(test_LambertTargeterIzzo "${SRCROOT}${MISSIONSEGMENTSDIR}")
#target_link_libraries(test_LambertTargeterIzzo tudat_mission_segments tudat_root_finders tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})
#
#add_executable(test_LambertTargeterGooding "${SRCROOT}${MISSIONSEGMENTSDIR}/tests/unitTestLambertTargeterGooding.cpp")
#setup_custom_test_program(test_LambertTargeterGooding "${SRCROOT}${MISSIONSEGMENTSDIR}")
#target_link_libraries(test_LambertTargeterGooding tudat_mission_segments tudat_root_finders tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})
#
#add_executable(test_LambertTargeter "${SRCROOT}${MISSIONSEGMENTSDIR}/tests/unitTestLambertTargeter.cpp")
#setup_custom_test_program(test_LambertTargeter "${SRCROOT}${MISSIONSEGMENTSDIR}")
#target_link_libraries(test_LambertTargeter tudat_mission_segments tudat_root_finders tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})
#
#add_executable(test_LambertRoutines "${SRCROOT}${MISSIONSEGMENTSDIR}/tests/unitTestLambertRoutines.cpp")
#setup_custom_test_program(test_LambertRoutines "${SRCROOT}${MISSIONSEGMENTSDIR}")
#target_link_libraries(test_LambertRoutines tudat_mission_segments tudat_root_finders tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})
#
#add_executable(test_ZeroRevolutionLambertTargeterIzzo "${SRCROOT}${MISSIONSEGMENTSDIR}/tests/unitTestZeroRevolutionLambertTargeterIzzo.cpp")
#setup_custom_test_program(test_ZeroRevolutionLambertTargeterIzzo "${SRCROOT}${MISSIONSEGMENTSDIR}")
#target_link_libraries(test_ZeroRevolutionLambertTargeterIzzo tudat_mission_segments tudat_root_finders tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})
#
#add_executable(test_MultiRevolutionLambertTargeterIzzo "${SRCROOT}${MISSIONSEGMENTSDIR}/tests/unitTestMultiRevolutionLambertTargeterIzzo.cpp")
#setup_custom_test_program(test_MultiRevolutionLambertTargeterIzzo "${SRCROOT}${MISSIONSEGMENTSDIR}")
#target_link_libraries(test_MultiRevolutionLambertTargeterIzzo tudat_mission_segments tudat_root_finders tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})
#
#add_executable(test_MathematicalShapeFunctions "${SRCROOT}${MISSIONSEGMENTSDIR}/tests/unitTestMathematicalShapeFunctions.cpp")
#setup_custom_test_program(test_MathematicalShapeFunctions "${SRCROOT}${MISSIONSEGMENTSDIR}")
#target_link_libraries(test_MathematicalShapeFunctions tudat_mission_segments tudat_basic_mathematics ${Boost_LIBRARIES})

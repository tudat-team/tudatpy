#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

TUDAT_ADD_TEST_CASE(EscapeAndCaptureRoutines PRIVATE_LINKS tudat_mission_segments tudat_gravitation tudat_basic_astrodynamics tudat_basic_mathematics tudat_root_finders tudat_input_output)

TUDAT_ADD_TEST_CASE(GravityAssistRoutines tudat_gravitation PRIVATE_LINKS tudat_mission_segments tudat_basic_astrodynamics tudat_basic_mathematics tudat_root_finders tudat_input_output)

TUDAT_ADD_TEST_CASE(LambertTargeterIzzo PRIVATE_LINKS tudat_mission_segments tudat_root_finders tudat_basic_astrodynamics tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(LambertTargeterGooding PRIVATE_LINKS tudat_mission_segments tudat_root_finders tudat_basic_astrodynamics tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(LambertTargeter PRIVATE_LINKS tudat_mission_segments tudat_root_finders tudat_basic_astrodynamics tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(LambertRoutines PRIVATE_LINKS tudat_mission_segments tudat_root_finders tudat_basic_astrodynamics tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(ZeroRevolutionLambertTargeterIzzo PRIVATE_LINKS tudat_mission_segments tudat_root_finders tudat_basic_astrodynamics tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(MultiRevolutionLambertTargeterIzzo PRIVATE_LINKS tudat_mission_segments tudat_root_finders tudat_basic_astrodynamics tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(MathematicalShapeFunctions PRIVATE_LINKS tudat_mission_segments tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(CaptureLeg PRIVATE_LINKS tudat_mission_segments tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(DepartureLegMga PRIVATE_LINKS tudat_mission_segments tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(MgaDsmPosition PRIVATE_LINKS tudat_mission_segments tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(MgaDsmVelocity PRIVATE_LINKS tudat_mission_segments tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(SwingbyLegMga PRIVATE_LINKS tudat_mission_segments tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(MgaTrajectory PRIVATE_LINKS tudat_mission_segments ${Tudat_PROPAGATION_LIBRARIES})


#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

TUDAT_ADD_TEST_CASE(CaptureLeg PRIVATE_LINKS tudat_trajectory_design tudat_mission_segments tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(DepartureLegMga PRIVATE_LINKS tudat_trajectory_design tudat_mission_segments tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(DepartureLegMga1DsmPosition PRIVATE_LINKS tudat_trajectory_design tudat_mission_segments tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(DepartureLegMga1DsmVelocity PRIVATE_LINKS tudat_trajectory_design tudat_mission_segments tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(SwingbyLegMga PRIVATE_LINKS tudat_trajectory_design tudat_mission_segments tudat_basic_mathematics)

#TUDAT_ADD_TEST_CASE(SwingbyLegMga1DsmPosition PRIVATE_LINKS tudat_trajectory_design tudat_mission_segments tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(SwingbyLegMga1DsmVelocity PRIVATE_LINKS tudat_trajectory_design tudat_mission_segments tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(Trajectory PRIVATE_LINKS tudat_trajectory_design tudat_mission_segments tudat_ephemerides tudat_basic_astrodynamics tudat_basic_mathematics)

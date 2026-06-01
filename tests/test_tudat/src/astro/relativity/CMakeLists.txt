#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

TUDAT_ADD_TEST_CASE(RelativisticAccelerationCorrection PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(EihAcceleration PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(ShapiroTimeDelay PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(RelativisticTimeConversion PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(SchwarzSchildMetric PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(SolarSystemMetric PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

# Temporarily removed: depends on local ISS CSV and writes to a user-specific path.
# TUDAT_ADD_TEST_CASE(IssProperTimeRateIau PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

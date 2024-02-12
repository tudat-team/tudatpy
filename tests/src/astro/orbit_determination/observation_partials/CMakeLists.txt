 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

TUDAT_ADD_TEST_CASE(OneWayRangePartials
        PRIVATE_LINKS
        tudat_test_support
        ${Tudat_ESTIMATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(NWayRangePartials
        PRIVATE_LINKS
        tudat_test_support
        ${Tudat_ESTIMATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(DifferencedOneWayRangePartials
        PRIVATE_LINKS
        tudat_test_support
        ${Tudat_ESTIMATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(DifferencedNWayRangePartials
        PRIVATE_LINKS
        tudat_test_support
        ${Tudat_ESTIMATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(OneWayDopplerPartials
        PRIVATE_LINKS
        tudat_test_support
        ${Tudat_ESTIMATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(TwoWayDopplerPartials
        PRIVATE_LINKS
        tudat_test_support
        ${Tudat_ESTIMATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(AngularPositionPartials
        PRIVATE_LINKS
        tudat_test_support
        ${Tudat_ESTIMATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(PositionPartials
        PRIVATE_LINKS
        tudat_test_support
        ${Tudat_ESTIMATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(EulerAngleObservationPartials
        PRIVATE_LINKS
        tudat_test_support
        ${Tudat_ESTIMATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(LightTimePartials
        PRIVATE_LINKS
        tudat_test_support
        ${Tudat_ESTIMATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(RotationMatrixPartials
        PRIVATE_LINKS
        ${Tudat_ESTIMATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(FullRotationPartials
        PRIVATE_LINKS
        ${Tudat_ESTIMATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(RelativeAngularPositionPartials
        PRIVATE_LINKS
        tudat_test_support
        ${Tudat_ESTIMATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(DsnNWayAveragedDopplerPartials
        PRIVATE_LINKS
        tudat_test_support
        ${Tudat_ESTIMATION_LIBRARIES}
        )
TUDAT_ADD_TEST_CASE(RelativePositionPartials
     PRIVATE_LINKS
     tudat_test_support
     ${Tudat_ESTIMATION_LIBRARIES}
     )
 TUDAT_ADD_TEST_CASE(TimeBiasPartials
         PRIVATE_LINKS
         tudat_test_support
         ${Tudat_ESTIMATION_LIBRARIES}
 )
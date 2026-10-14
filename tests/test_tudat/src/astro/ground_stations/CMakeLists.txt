#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

TUDAT_ADD_TEST_CASE(GroundStationState
    PRIVATE_LINKS
    ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(StationHistoryTrajectory
    PRIVATE_LINKS
    ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(StationHistoryFrameTransformation
    PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(StationReferencePoint
    PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(PointingAnglesCalculator
    PRIVATE_LINKS
    ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(DefaultGroundStationSettings
    PRIVATE_LINKS
    ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(Iers2010DeformationModel
        PRIVATE_LINKS
        ${Tudat_ESTIMATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(OceanTideDisplacement
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
)

TUDAT_ADD_TEST_CASE(PoleTideDisplacement
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
)

 TUDAT_ADD_TEST_CASE(TransmittingFrequencies
     PRIVATE_LINKS
     ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(GroundStationSunVisibility
    PRIVATE_LINKS
    ${Tudat_ESTIMATION_LIBRARIES})
TUDAT_ADD_TEST_CASE(OceanPoleTideDisplacement PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(BodyDeformationVelocity
    PRIVATE_LINKS
    ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(EarthDeformationTimeRouting
        PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(AtmosphericTideLoading PRIVATE_LINKS tudat_ground_stations tudat_environment_data)

TUDAT_ADD_TEST_CASE(AtmosphericTideDeformation PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES} tudat_environment_data)

TUDAT_ADD_TEST_CASE(CompleteEarthStationVelocity
    PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES} tudat_environment_data)

TUDAT_ADD_TEST_CASE(RelativisticStationPosition PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(StationLoadingDeformation PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

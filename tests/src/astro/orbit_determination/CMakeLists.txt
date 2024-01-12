 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

#set(ACCELERATIONPARTIALSDIR "${ORBITDETERMINATIONDIR}/AccelerationPartials")
#set(TORQUEPARTIALSDIR "${ORBITDETERMINATIONDIR}/RotationalDynamicsPartials")
#set(ESTIMATABLEPARAMETERSDIR "${ORBITDETERMINATIONDIR}/EstimatableParameters")
#set(OBSERVATIONPARTIALSDIR "${ORBITDETERMINATIONDIR}/ObservationPartials")
#set(LIGHTTIMECORRECTIONPARTIALSDIR "${ORBITDETERMINATIONDIR}/LightTimeCorrectionPartials")

# Add subdirectories.
add_subdirectory(acceleration_partials)
add_subdirectory(rotational_dynamics_partials)
add_subdirectory(observation_partials)

TUDAT_ADD_TEST_CASE(EstimationInput
    PRIVATE_LINKS
    ${Tudat_ESTIMATION_LIBRARIES}
    )

TUDAT_ADD_TEST_CASE(EstimationFromIdealDataDoubleDouble
    PRIVATE_LINKS
    ${Tudat_ESTIMATION_LIBRARIES}
    )

TUDAT_ADD_TEST_CASE(TidalPropertyEstimation
    PRIVATE_LINKS
    ${Tudat_ESTIMATION_LIBRARIES}
    )

TUDAT_ADD_TEST_CASE(HybridArcStateEstimation
    PRIVATE_LINKS
    ${Tudat_ESTIMATION_LIBRARIES}
    )

# TUDAT_ADD_TEST_CASE(HybridArcMultiBodyEstimation
#         PRIVATE_LINKS
#         ${Tudat_ESTIMATION_LIBRARIES}
#         )

TUDAT_ADD_TEST_CASE(MultiArcStateEstimation
    PRIVATE_LINKS
    ${Tudat_ESTIMATION_LIBRARIES}
    )

TUDAT_ADD_TEST_CASE(ParameterInfluenceDetermination
    PRIVATE_LINKS
    ${Tudat_ESTIMATION_LIBRARIES}
    )

TUDAT_ADD_TEST_CASE(BiasEstimation
    PRIVATE_LINKS
    ${Tudat_ESTIMATION_LIBRARIES}
    )

TUDAT_ADD_TEST_CASE(TimeBiasEstimation
     PRIVATE_LINKS
     ${Tudat_ESTIMATION_LIBRARIES}
     )

TUDAT_ADD_TEST_CASE(RotationalStateEstimation
    PRIVATE_LINKS
    ${Tudat_ESTIMATION_LIBRARIES}
    )

TUDAT_ADD_TEST_CASE(ArcwiseEnvironmentParameters
    PRIVATE_LINKS
    ${Tudat_ESTIMATION_LIBRARIES}
    )

TUDAT_ADD_TEST_CASE(DesaturationDeltaVsEstimation
    PRIVATE_LINKS
    ${Tudat_ESTIMATION_LIBRARIES}
    )

TUDAT_ADD_TEST_CASE(FullPlanetaryRotationModelEstimation
    PRIVATE_LINKS
    ${Tudat_ESTIMATION_LIBRARIES}
    )


TUDAT_ADD_TEST_CASE(NonSequentialStateEstimation
     PRIVATE_LINKS
     ${Tudat_ESTIMATION_LIBRARIES}
     )

 TUDAT_ADD_TEST_CASE(EihEstimation
         PRIVATE_LINKS
         ${Tudat_ESTIMATION_LIBRARIES}
         )
# Example application: estimation of Doppler observations from ODF/DSN data
#TUDAT_ADD_TEST_CASE(EstimationFromDsnData
#    PRIVATE_LINKS
#    ${Tudat_ESTIMATION_LIBRARIES}
#    )

# Example application: estimation of position observables from SPICE data
#TUDAT_ADD_TEST_CASE(EstimationFromSpiceData
#    PRIVATE_LINKS
#    ${Tudat_ESTIMATION_LIBRARIES}
#    )

 TUDAT_ADD_TEST_CASE(ConsiderParameters
         PRIVATE_LINKS
         ${Tudat_ESTIMATION_LIBRARIES}
         )

 TUDAT_ADD_TEST_CASE(FitToSpice
        PRIVATE_LINKS
        ${Tudat_ESTIMATION_LIBRARIES}
 )

 TUDAT_ADD_TEST_CASE(CustomAccelerationVariationalEquations
         PRIVATE_LINKS
         ${Tudat_ESTIMATION_LIBRARIES}
 )


 TUDAT_ADD_TEST_CASE(VehicleReferencePointEstimation
         PRIVATE_LINKS
         ${Tudat_ESTIMATION_LIBRARIES}
 )

 TUDAT_ADD_TEST_CASE(RadiationPressureEstimation
         PRIVATE_LINKS
         ${Tudat_ESTIMATION_LIBRARIES}
 )

if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )


    TUDAT_ADD_TEST_CASE(EstimationFromIdealDataTimeDouble
            PRIVATE_LINKS
            ${Tudat_ESTIMATION_LIBRARIES}
    )

    TUDAT_ADD_TEST_CASE(EstimationFromIdealDataDoubleLongDouble
        PRIVATE_LINKS
        ${Tudat_ESTIMATION_LIBRARIES}
        )

    TUDAT_ADD_TEST_CASE(EstimationFromIdealDataTimeLongDouble
        PRIVATE_LINKS
        ${Tudat_ESTIMATION_LIBRARIES}
        )

endif( )



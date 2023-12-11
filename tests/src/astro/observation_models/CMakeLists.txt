#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.

TUDAT_ADD_TEST_CASE(LightTimeSolution
        PRIVATE_LINKS
        tudat_observation_models
        tudat_ephemerides
        tudat_spice_interface
        )


TUDAT_ADD_TEST_CASE(AngularPositionModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})


TUDAT_ADD_TEST_CASE(OneWayRangeModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})


TUDAT_ADD_TEST_CASE(OneWayDifferencedRangeModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})


TUDAT_ADD_TEST_CASE(DopplerModels PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})


TUDAT_ADD_TEST_CASE(PositionObservationModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})


TUDAT_ADD_TEST_CASE(EulerAngleObservationModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})


TUDAT_ADD_TEST_CASE(VelocityObservationModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})


TUDAT_ADD_TEST_CASE(SimulatedObservationNoise PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})


TUDAT_ADD_TEST_CASE(NWayRangeObservationModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})


TUDAT_ADD_TEST_CASE(NWayRangeRateObservationModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})


TUDAT_ADD_TEST_CASE(RelativeAngularPositionModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(RelativePositionModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})


TUDAT_ADD_TEST_CASE(ObservationViabilityCalculators PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})


TUDAT_ADD_TEST_CASE(PerArcObservationSimulation PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})


TUDAT_ADD_TEST_CASE(TimeBias PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(AtmosphereCorrection PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(SolarCoronaCorrection PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(ObservationRoundingErrors PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(ObservationDependentVariables PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})


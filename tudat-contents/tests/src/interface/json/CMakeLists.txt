#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#


TUDAT_ADD_TEST_CASE(Acceleration PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(Aerodynamics PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(Atmosphere PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(Body PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(Deserialization PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(Ephemeris PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(Export PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(GravityField PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(GravityFieldVariation PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(GroundStation PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(Integrator PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(Interpolation PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(MassRateModel PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(Propagator PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(RadiationPressure PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(ReferenceFrames PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(RotationModel PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(ShapeModel PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(SimulationSingleSatellite PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(SimulationSinglePerturbedSatellite PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(SimulationInnerSolarSystem PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(SimulationGalileoConstellation PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(SimulationThrustAlongVelocityVector PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(SimulationThrustAccelerationFromFile PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(Spice PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(State PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(Termination PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(Thrust PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(Torque PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(Variable PRIVATE_LINKS tudat_json_interface_library ${Tudat_PROPAGATION_LIBRARIES})

if (BUILD_WITH_ESTIMATION_TOOLS)
    TUDAT_ADD_TEST_CASE(Observation PRIVATE_LINKS tudat_json_interface_library ${Tudat_ESTIMATION_LIBRARIES})

    TUDAT_ADD_TEST_CASE(Parameter PRIVATE_LINKS tudat_json_interface_library ${Tudat_ESTIMATION_LIBRARIES})

    TUDAT_ADD_TEST_CASE(SimulationSingleSatelliteVariational PRIVATE_LINKS tudat_json_interface_library ${Tudat_ESTIMATION_LIBRARIES})

    TUDAT_ADD_TEST_CASE(SimulationSingleSatelliteEstimation PRIVATE_LINKS tudat_json_interface_library ${Tudat_ESTIMATION_LIBRARIES})
endif ()



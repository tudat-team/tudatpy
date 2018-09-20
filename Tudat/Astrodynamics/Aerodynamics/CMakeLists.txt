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
set(AERODYNAMICS_SOURCES
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicAcceleration.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicCoefficientGenerator.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicTorque.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicForce.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamics.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/customConstantTemperatureAtmosphere.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/exponentialAtmosphere.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/hypersonicLocalInclinationAnalysis.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/tabulatedAtmosphere.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/flightConditions.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/trimOrientation.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/equilibriumWallTemperature.cpp"
)

# Set the header files.
set(AERODYNAMICS_HEADERS
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicAcceleration.h"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicTorque.h"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicCoefficientGenerator.h"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicCoefficientInterface.h"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicAcceleration.h"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicForce.h"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamics.h"
  "${SRCROOT}${AERODYNAMICSDIR}/atmosphereModel.h"
  "${SRCROOT}${AERODYNAMICSDIR}/customConstantTemperatureAtmosphere.h"
  "${SRCROOT}${AERODYNAMICSDIR}/exponentialAtmosphere.h"
  "${SRCROOT}${AERODYNAMICSDIR}/hypersonicLocalInclinationAnalysis.h"
  "${SRCROOT}${AERODYNAMICSDIR}/tabulatedAtmosphere.h"
  "${SRCROOT}${AERODYNAMICSDIR}/standardAtmosphere.h"
  "${SRCROOT}${AERODYNAMICSDIR}/customAerodynamicCoefficientInterface.h"
  "${SRCROOT}${AERODYNAMICSDIR}/controlSurfaceAerodynamicCoefficientInterface.h"
  "${SRCROOT}${AERODYNAMICSDIR}/flightConditions.h"
  "${SRCROOT}${AERODYNAMICSDIR}/UnitTests/testApolloCapsuleCoefficients.h"
  "${SRCROOT}${AERODYNAMICSDIR}/flightConditions.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicGuidance.h"
  "${SRCROOT}${AERODYNAMICSDIR}/equilibriumWallTemperature.h"
  "${SRCROOT}${AERODYNAMICSDIR}/windModel.h"
)

if(USE_NRLMSISE00)
  set(AERODYNAMICS_SOURCES "${AERODYNAMICS_SOURCES}"
    "${SRCROOT}${AERODYNAMICSDIR}/nrlmsise00Atmosphere.cpp"
    "${SRCROOT}${AERODYNAMICSDIR}/nrlmsise00InputFunctions.cpp")
  set(AERODYNAMICS_HEADERS "${AERODYNAMICS_HEADERS}"
    "${SRCROOT}${AERODYNAMICSDIR}/nrlmsise00Atmosphere.h"
    "${SRCROOT}${AERODYNAMICSDIR}/nrlmsise00InputFunctions.h")
endif( )

# Add static libraries.
add_library(tudat_aerodynamics STATIC ${AERODYNAMICS_SOURCES} ${AERODYNAMICS_HEADERS})
setup_tudat_library_target(tudat_aerodynamics "${SRCROOT}{AERODYNAMICSDIR}")

# Add unit tests.
add_executable(test_AerodynamicMomentAndAerodynamicForce "${SRCROOT}${AERODYNAMICSDIR}/UnitTests/unitTestAerodynamicMomentAndAerodynamicForce.cpp")
setup_custom_test_program(test_AerodynamicMomentAndAerodynamicForce "${SRCROOT}${AERODYNAMICSDIR}")
target_link_libraries(test_AerodynamicMomentAndAerodynamicForce ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_AerodynamicsNamespace "${SRCROOT}${AERODYNAMICSDIR}/UnitTests/unitTestAerodynamicsNamespace.cpp")
setup_custom_test_program(test_AerodynamicsNamespace "${SRCROOT}${AERODYNAMICSDIR}")
target_link_libraries(test_AerodynamicsNamespace tudat_aerodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_AerodynamicCoefficientGenerator "${SRCROOT}${AERODYNAMICSDIR}/UnitTests/unitTestCoefficientGenerator.cpp")
setup_custom_test_program(test_AerodynamicCoefficientGenerator "${SRCROOT}${AERODYNAMICSDIR}")
target_link_libraries(test_AerodynamicCoefficientGenerator tudat_aerodynamics tudat_geometric_shapes tudat_interpolators tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_ExponentialAtmosphere "${SRCROOT}${AERODYNAMICSDIR}/UnitTests/unitTestExponentialAtmosphere.cpp")
setup_custom_test_program(test_ExponentialAtmosphere "${SRCROOT}${AERODYNAMICSDIR}")
target_link_libraries(test_ExponentialAtmosphere tudat_aerodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_CustomConstantTemperatureAtmosphere "${SRCROOT}${AERODYNAMICSDIR}/UnitTests/unitTestCustomConstantTemperatureAtmosphere.cpp")
setup_custom_test_program(test_CustomConstantTemperatureAtmosphere "${SRCROOT}${AERODYNAMICSDIR}")
target_link_libraries(test_CustomConstantTemperatureAtmosphere tudat_aerodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_TabulatedAtmosphere "${SRCROOT}${AERODYNAMICSDIR}/UnitTests/unitTestTabulatedAtmosphere.cpp")
setup_custom_test_program(test_TabulatedAtmosphere "${SRCROOT}${AERODYNAMICSDIR}")
target_link_libraries(test_TabulatedAtmosphere tudat_aerodynamics tudat_input_output tudat_interpolators
    tudat_basic_mathematics tudat_basic_astrodynamics ${Boost_LIBRARIES})

add_executable(test_TabulatedAerodynamicCoefficients "${SRCROOT}${AERODYNAMICSDIR}/UnitTests/unitTestTabulatedAerodynamicCoefficients.cpp")
setup_custom_test_program(test_TabulatedAerodynamicCoefficients "${SRCROOT}${AERODYNAMICSDIR}")
target_link_libraries(test_TabulatedAerodynamicCoefficients ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_HeatTransfer "${SRCROOT}${AERODYNAMICSDIR}/UnitTests/unitTestHeatTransfer.cpp")
setup_custom_test_program(test_HeatTransfer "${SRCROOT}${AERODYNAMICSDIR}")
target_link_libraries(test_HeatTransfer tudat_aerodynamics tudat_root_finders tudat_basic_mathematics
    tudat_basic_astrodynamics ${Boost_LIBRARIES})

add_executable(test_ControlSurfaceIncrements "${SRCROOT}${AERODYNAMICSDIR}/UnitTests/unitTestControlSurfaceIncrements.cpp")
setup_custom_test_program(test_ControlSurfaceIncrements "${SRCROOT}${AERODYNAMICSDIR}")
target_link_libraries(test_ControlSurfaceIncrements ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_AerodynamicCoefficientsFromFile "${SRCROOT}${AERODYNAMICSDIR}/UnitTests/unitTestAerodynamicCoefficientsFromFile.cpp")
setup_custom_test_program(test_AerodynamicCoefficientsFromFile "${SRCROOT}${AERODYNAMICSDIR}")
target_link_libraries(test_AerodynamicCoefficientsFromFile ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_WindModel "${SRCROOT}${AERODYNAMICSDIR}/UnitTests/unitTestWindModel.cpp")
setup_custom_test_program(test_WindModel "${SRCROOT}${AERODYNAMICSDIR}")
target_link_libraries(test_WindModel ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

if(USE_NRLMSISE00)
    add_executable(test_NRLMSISE00Atmosphere "${SRCROOT}${AERODYNAMICSDIR}/UnitTests/unitTestNRLMSISE00Atmosphere.cpp")
    setup_custom_test_program(test_NRLMSISE00Atmosphere "${SRCROOT}${AERODYNAMICSDIR}")
    target_link_libraries(test_NRLMSISE00Atmosphere tudat_aerodynamics tudat_interpolators tudat_basic_mathematics nrlmsise00 tudat_input_output tudat_basic_astrodynamics ${Boost_LIBRARIES})
endif( )

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
set(electromagnetism_HEADERS
        "cannonBallRadiationPressureAcceleration.h"
        "cannonBallRadiationPressureForce.h"
        "lorentzStaticMagneticForce.h"
        "lorentzStaticMagneticAcceleration.h"
        "radiationPressureInterface.h"
        "basicElectroMagnetism.h"
        "panelledRadiationPressure.h"
        "solarSailAcceleration.h"
        "solarSailForce.h"
        )

# Set the header files.
set(electromagnetism_SOURCES
        "cannonBallRadiationPressureAcceleration.cpp"
        "cannonBallRadiationPressureForce.cpp"
        "lorentzStaticMagneticForce.cpp"
        "lorentzStaticMagneticAcceleration.cpp"
        "radiationPressureInterface.cpp"
        "panelledRadiationPressure.cpp"
        "solarSailAcceleration.cpp"
        "solarSailForce.cpp"
        )

# Add tests.
TUDAT_ADD_LIBRARY("electromagnetism"
        "${electromagnetism_SOURCES}"
        "${electromagnetism_HEADERS}"
        PRIVATE_LINKS "${Boost_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )

TUDAT_ADD_TEST_CASE(CannonBallRadiationPressureAccelerationAndForce
        PRIVATE_LINKS
        tudat_electromagnetism
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(LorentzStaticMagneticAccelerationAndForce
        PRIVATE_LINKS
        tudat_electromagnetism
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(RadiationPressureInterface
        PRIVATE_LINKS
        tudat_electromagnetism
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(PanelledRadiationPressure
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
        tudat_electromagnetism
        tudat_spice_interface
        ${CSpice_LIBRARIES}
        ${Boost_LIBRARIES}
        ${Sofa_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(SolarSailAccelerationAndForce
        PRIVATE_LINKS
        tudat_electromagnetism
        tudat_basic_astrodynamics
        )
# Add unit tests.
#add_executable(test_CannonBallRadiationPressureAccelerationAndForce "${SRCROOT}${ELECTROMAGNETISMDIR}/tests/unitTestCannonBallRadiationPressureAccelerationAndForce.cpp")
#setup_custom_test_program(test_CannonBallRadiationPressureAccelerationAndForce "${SRCROOT}${ELECTROMAGNETISMDIR}")
#target_link_libraries(test_CannonBallRadiationPressureAccelerationAndForce tudat_electromagnetism tudat_basic_astrodynamics ${Boost_LIBRARIES})
#
#add_executable(test_LorentzStaticMagneticAccelerationAndForce "${SRCROOT}${ELECTROMAGNETISMDIR}/tests/unitTestLorentzStaticMagneticAccelerationAndForce.cpp")
#setup_custom_test_program(test_LorentzStaticMagneticAccelerationAndForce "${SRCROOT}${ELECTROMAGNETISMDIR}")
#target_link_libraries(test_LorentzStaticMagneticAccelerationAndForce tudat_electromagnetism tudat_basic_astrodynamics ${Boost_LIBRARIES})
#
#add_executable(test_RadiationPressureInterface "${SRCROOT}${ELECTROMAGNETISMDIR}/tests/unitTestRadiationPressureInterface.cpp")
#setup_custom_test_program(test_RadiationPressureInterface "${SRCROOT}${ELECTROMAGNETISMDIR}")
#target_link_libraries(test_RadiationPressureInterface tudat_electromagnetism tudat_basic_astrodynamics ${Boost_LIBRARIES})
#
#add_executable(test_PanelledRadiationPressure "${SRCROOT}${ELECTROMAGNETISMDIR}/tests/unitTestPanelledRadiationPressure.cpp")
#setup_custom_test_program(test_PanelledRadiationPressure "${ELECTROMAGNETISMDIR}")
#target_link_libraries(test_PanelledRadiationPressure ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_SolarSailAccelerationAndForce "${SRCROOT}${ELECTROMAGNETISMDIR}/tests/unitTestSolarSailAccelerationAndForce.cpp")
#setup_custom_test_program(test_SolarSailAccelerationAndForce "${SRCROOT}${ELECTROMAGNETISMDIR}")
#target_link_libraries(test_SolarSailAccelerationAndForce tudat_electromagnetism tudat_basic_astrodynamics ${Boost_LIBRARIES})

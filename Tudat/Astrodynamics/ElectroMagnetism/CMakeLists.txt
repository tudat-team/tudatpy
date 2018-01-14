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
set(ELECTROMAGNETISM_SOURCES
  "${SRCROOT}${ELECTROMAGNETISMDIR}/cannonBallRadiationPressureAcceleration.h"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/cannonBallRadiationPressureForce.h"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/lorentzStaticMagneticForce.h"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/lorentzStaticMagneticAcceleration.h"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/radiationPressureInterface.h"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/basicElectroMagnetism.h"
)

# Set the header files.
set(ELECTROMAGNETISM_HEADERS
  "${SRCROOT}${ELECTROMAGNETISMDIR}/cannonBallRadiationPressureAcceleration.cpp"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/cannonBallRadiationPressureForce.cpp"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/lorentzStaticMagneticForce.cpp"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/lorentzStaticMagneticAcceleration.cpp"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/radiationPressureInterface.cpp"
)

# Add static libraries.
add_library(tudat_electro_magnetism STATIC ${ELECTROMAGNETISM_SOURCES} ${ELECTROMAGNETISM_HEADERS})
setup_tudat_library_target(tudat_electro_magnetism "${SRCROOT}${ELECTROMAGNETISMDIR}")

# Add unit tests.
add_executable(test_CannonBallRadiationPressureAccelerationAndForce "${SRCROOT}${ELECTROMAGNETISMDIR}/UnitTests/unitTestCannonBallRadiationPressureAccelerationAndForce.cpp")
setup_custom_test_program(test_CannonBallRadiationPressureAccelerationAndForce "${SRCROOT}${ELECTROMAGNETISMDIR}")
target_link_libraries(test_CannonBallRadiationPressureAccelerationAndForce tudat_electro_magnetism ${Boost_LIBRARIES})

add_executable(test_LorentzStaticMagneticAccelerationAndForce "${SRCROOT}${ELECTROMAGNETISMDIR}/UnitTests/unitTestLorentzStaticMagneticAccelerationAndForce.cpp")
setup_custom_test_program(test_LorentzStaticMagneticAccelerationAndForce "${SRCROOT}${ELECTROMAGNETISMDIR}")
target_link_libraries(test_LorentzStaticMagneticAccelerationAndForce tudat_electro_magnetism ${Boost_LIBRARIES})

add_executable(test_RadiationPressureInterface "${SRCROOT}${ELECTROMAGNETISMDIR}/UnitTests/unitTestRadiationPressureInterface.cpp")
setup_custom_test_program(test_RadiationPressureInterface "${SRCROOT}${ELECTROMAGNETISMDIR}")
target_link_libraries(test_RadiationPressureInterface tudat_electro_magnetism tudat_basic_astrodynamics ${Boost_LIBRARIES})

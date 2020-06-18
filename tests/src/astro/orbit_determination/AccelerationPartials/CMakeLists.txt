#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.

if(TUDAT_BUILD_WITH_SPICE_INTERFACE)

# Add unit tests
add_executable(test_AccelerationPartials "${SRCROOT}${ACCELERATIONPARTIALSDIR}/tests/unitTestAccelerationPartials.cpp")
setup_custom_test_program(test_AccelerationPartials "${SRCROOT}${ACCELERATIONPARTIALSDIR}")
target_link_libraries(test_AccelerationPartials ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

# Add unit tests
add_executable(test_SphericalHarmonicPartials "${SRCROOT}${ACCELERATIONPARTIALSDIR}/tests/unitTestSphericalHarmonicPartials.cpp")
setup_custom_test_program(test_SphericalHarmonicPartials "${SRCROOT}${ACCELERATIONPARTIALSDIR}")
target_link_libraries(test_SphericalHarmonicPartials ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_MutualSphericalHarmonicPartials "${SRCROOT}${ACCELERATIONPARTIALSDIR}/tests/unitTestMutualSphericalHarmonicPartials.cpp")
setup_custom_test_program(test_MutualSphericalHarmonicPartials "${SRCROOT}${ACCELERATIONPARTIALSDIR}")
target_link_libraries(test_MutualSphericalHarmonicPartials ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

endif( )



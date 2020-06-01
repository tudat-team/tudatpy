#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# http://www.iausofa.org/current_C.html#Downloads

#Set the source files.
set(sofa_interface_SOURCES
        "sofaTimeConversions.cpp"
        "fundamentalArguments.cpp"
        "earthOrientation.cpp"
        )

# Set the header files.
set(sofa_interface_HEADERS
        "sofaTimeConversions.h"
        "fundamentalArguments.h"
        "earthOrientation.h"
        )

# Add static libraries.
TUDAT_ADD_LIBRARY("sofa_interface"
        "${sofa_interface_SOURCES}"
        "${sofa_interface_HEADERS}"
        PUBLIC_LINKS
        tudat_basic_astrodynamics
        PRIVATE_LINKS "${Boost_LIBRARIES}" "${Sofa_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}" "${Sofa_INCLUDE_DIRS}"
        )

TUDAT_ADD_TEST_CASE(SofaEarthOrientationInterface
        PRIVATE_LINKS
        tudat_sofa_interface
        tudat_basic_mathematics
        "${Sofa_LIBRARIES}"
        )

TUDAT_ADD_TEST_CASE(SofaFundamentalArguments
        PRIVATE_LINKS
        tudat_sofa_interface
        tudat_basic_mathematics
        tudat_basic_astrodynamics
        "${Sofa_LIBRARIES}"

        )

TUDAT_ADD_TEST_CASE(SofaTimeConversions
        PRIVATE_LINKS
        tudat_sofa_interface
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        "${Sofa_LIBRARIES}"
        )


#add_executable(test_SofaEarthOrientation "${SRCROOT}${EXTERNALDIR}/SofaInterface/tests/unitTestSofaEarthOrientationInterface.cpp")
#setup_custom_test_program(test_SofaEarthOrientation "${SRCROOT}${EXTERNALDIR}/SofaInterface")
#target_link_libraries(test_SofaEarthOrientation tudat_sofa_interface tudat_basic_mathematics sofa ${Boost_LIBRARIES})
#
#add_executable(test_SofaFundamentalArguments "${SRCROOT}${EXTERNALDIR}/SofaInterface/tests/unitTestFundamentalArguments.cpp")
#setup_custom_test_program(test_SofaFundamentalArguments "${SRCROOT}${EXTERNALDIR}/SofaInterface")
#target_link_libraries(test_SofaFundamentalArguments tudat_sofa_interface tudat_basic_astrodynamics sofa ${Boost_LIBRARIES})
#
#add_executable(test_SofaTimeConversions "${SRCROOT}${EXTERNALDIR}/SofaInterface/tests/unitTestSofaTimeConversions.cpp")
#setup_custom_test_program(test_SofaTimeConversions "${SRCROOT}${EXTERNALDIR}/SofaInterface")
#target_link_libraries(test_SofaTimeConversions tudat_basic_astrodynamics tudat_basic_mathematics tudat_sofa_interface sofa ${Boost_LIBRARIES})

# End of file.

#    Copyright (c) 2010-2018, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.

#Set the source files.
set(SOFAINTERFACE_SOURCES
  "${SRCROOT}${EXTERNALDIR}/SofaInterface/sofaTimeConversions.h"
  "${SRCROOT}${EXTERNALDIR}/SofaInterface/fundamentalArguments.h"
  "${SRCROOT}${EXTERNALDIR}/SofaInterface/earthOrientation.h"

)

# Set the header files.
set(SOFAINTERFACE_HEADERS
"${SRCROOT}${EXTERNALDIR}/SofaInterface/sofaTimeConversions.cpp"
"${SRCROOT}${EXTERNALDIR}/SofaInterface/fundamentalArguments.cpp"
  "${SRCROOT}${EXTERNALDIR}/SofaInterface/earthOrientation.cpp"
)

# Add static libraries.
add_library(tudat_sofa_interface STATIC ${SOFAINTERFACE_SOURCES} ${SOFAINTERFACE_HEADERS})
setup_tudat_library_target(tudat_sofa_interface "${SRCROOT}${SOFAINTERFACEDIR}")

add_executable(test_SofaEarthOrientation "${SRCROOT}${EXTERNALDIR}/SofaInterface/UnitTests/unitTestSofaEarthOrientationInterface.cpp")
setup_custom_test_program(test_SofaEarthOrientation "${SRCROOT}${EXTERNALDIR}/SofaInterface")
target_link_libraries(test_SofaEarthOrientation tudat_sofa_interface tudat_basic_mathematics sofa ${Boost_LIBRARIES})

add_executable(test_SofaFundamentalArguments "${SRCROOT}${EXTERNALDIR}/SofaInterface/UnitTests/unitTestFundamentalArguments.cpp")
setup_custom_test_program(test_SofaFundamentalArguments "${SRCROOT}${EXTERNALDIR}/SofaInterface")
target_link_libraries(test_SofaFundamentalArguments tudat_sofa_interface tudat_basic_astrodynamics sofa ${Boost_LIBRARIES})

add_executable(test_SofaTimeConversions "${SRCROOT}${EXTERNALDIR}/SofaInterface/UnitTests/unitTestSofaTimeConversions.cpp")
setup_custom_test_program(test_SofaTimeConversions "${SRCROOT}${EXTERNALDIR}/SofaInterface")
target_link_libraries(test_SofaTimeConversions tudat_basic_astrodynamics tudat_basic_mathematics tudat_sofa_interface sofa ${Boost_LIBRARIES})

# End of file.

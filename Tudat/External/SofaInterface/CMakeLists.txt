 #     Copyright (c) 2010-2012 Delft University of Technology.
 #
 #    This software is protected by national and international copyright.
 #    Any unauthorized use, reproduction or modification is unlawful and
 #    will be prosecuted. Commercial and non-private application of the
 #    software in any form is strictly prohibited unless otherwise granted
 #    by the authors.
 #
 #    The code is provided without any warranty; without even the implied
 #    warranty of merchantibility or fitness for a particular purpose.
 #
 #    Changelog
 #      YYMMDD    Author            Comment
 #      110820    S.M. Persson      File created.
 #      111025    K. Kumar          Adapted file to work with Revision 194.
 #      111026    K. Kumar          Adapted file so all headers show in project tree in Qt Creator.
 #/

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


#add_executable(test_UtcTdbConversion "${SRCROOT}${EXTERNALDIR}/SofaInterface/UnitTests/unitTestUtcTdbConversion.cpp")
#setup_custom_test_program(test_UtcTdbConversion "${SRCROOT}${EXTERNALDIR}/SofaInterface")
#target_link_libraries(test_UtcTdbConversion ilr_input_output tudat_input_output ilr_basic_astrodynamics tudat_basic_astrodynamics tudat_sofa_interface sofa gomp ${Boost_LIBRARIES})


# End of file.

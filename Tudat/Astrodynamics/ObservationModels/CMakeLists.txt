 #    Copyright (c) 2010-2012 Delft University of Technology.
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
 #      130226    D. Dirkx          File created.
 #

# Set the source files.
set(OBSERVATION_MODELS_SOURCES
  "${SRCROOT}${OBSERVATIONMODELSDIR}/lightTimeSolution.cpp"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/linkTypeDefs.cpp"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/observableTypes.cpp"
)

# Set the header files.
set(OBSERVATION_MODELS_HEADERS
  "${SRCROOT}${OBSERVATIONMODELSDIR}/angularPositionObservationModel.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/lightTimeSolution.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/linkTypeDefs.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/observableTypes.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/observationModel.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/oneWayRangeObservationModel.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/nWayRangeObservationModel.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/UnitTests/testLightTimeCorrections.h"
)


# Add static libraries.
add_library(tudat_observation_models STATIC ${OBSERVATION_MODELS_SOURCES} ${OBSERVATION_MODELS_HEADERS})
setup_tudat_library_target(tudat_observation_models "${SRCROOT}${OBSERVATIONMODELSDIR}")

if(USE_CSPICE)

    add_executable(test_LightTime "${SRCROOT}${OBSERVATIONMODELSDIR}/UnitTests/unitTestLightTimeSolution.cpp")
    setup_custom_test_program(test_LightTime "${SRCROOT}${OBSERVATIONMODELSDIR}")
    target_link_libraries(test_LightTime tudat_observation_models tudat_ephemerides tudat_spice_interface ${SPICE_LIBRARIES} ${Boost_LIBRARIES})

endif( )

    add_executable(test_AngularPositionModel "${SRCROOT}${OBSERVATIONMODELSDIR}/UnitTests/unitTestAngularPositionModel.cpp")
    setup_custom_test_program(test_AngularPositionModel "${SRCROOT}${OBSERVATIONMODELSDIR}")
    target_link_libraries(test_AngularPositionModel tudat_observation_models tudat_ephemerides tudat_spice_interface ${SPICE_LIBRARIES} ${Boost_LIBRARIES})

    add_executable(test_NWayRangePositionModel "${SRCROOT}${OBSERVATIONMODELSDIR}/UnitTests/unitTestNWayObservationModels.cpp")
    setup_custom_test_program(test_NWayRangePositionModel "${SRCROOT}${OBSERVATIONMODELSDIR}")
    target_link_libraries(test_NWayRangePositionModel tudat_observation_models tudat_ephemerides tudat_spice_interface ${SPICE_LIBRARIES} ${Boost_LIBRARIES})

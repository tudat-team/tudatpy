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


set(OBSERVABLECORRECTIONSDIR "${OBSERVATIONMODELSDIR}/ObservableCorrections")


# Set the source files.
set(OBSERVATION_MODELS_SOURCES
  "${SRCROOT}${OBSERVATIONMODELSDIR}/lightTimeSolution.cpp"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/observableTypes.cpp"
  "${SRCROOT}${OBSERVABLECORRECTIONSDIR}/createLightTimeCorrection.cpp"
  "${SRCROOT}${OBSERVABLECORRECTIONSDIR}/firstOrderRelativisticLightTimeCorrection.cpp"
)

# Set the header files.
set(OBSERVATION_MODELS_HEADERS
  "${SRCROOT}${OBSERVATIONMODELSDIR}/angularPositionObservationModel.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/lightTimeSolution.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/linkTypeDefs.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/observableTypes.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/observationModel.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/observationBias.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/oneWayRangeObservationModel.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/positionObservationModel.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/createObservationModel.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/createLightTimeCalculator.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/UnitTests/testLightTimeCorrections.h"
  "${SRCROOT}${OBSERVABLECORRECTIONSDIR}/createLightTimeCorrection.h"
  "${SRCROOT}${OBSERVABLECORRECTIONSDIR}/lightTimeCorrection.h"
  "${SRCROOT}${OBSERVABLECORRECTIONSDIR}/firstOrderRelativisticLightTimeCorrection.h"
)


# Add static libraries.
add_library(tudat_observation_models STATIC ${OBSERVATION_MODELS_SOURCES} ${OBSERVATION_MODELS_HEADERS})
setup_tudat_library_target(tudat_observation_models "${SRCROOT}${OBSERVATIONMODELSDIR}")

if(USE_CSPICE)

    add_executable(test_LightTime "${SRCROOT}${OBSERVATIONMODELSDIR}/UnitTests/unitTestLightTimeSolution.cpp")
    setup_custom_test_program(test_LightTime "${SRCROOT}${OBSERVATIONMODELSDIR}")
    target_link_libraries(test_LightTime tudat_observation_models tudat_ephemerides tudat_basic_astrodynamics tudat_spice_interface ${SPICE_LIBRARIES} ${Boost_LIBRARIES})

endif( )

    add_executable(test_AngularPositionModel "${SRCROOT}${OBSERVATIONMODELSDIR}/UnitTests/unitTestAngularPositionModel.cpp")
    setup_custom_test_program(test_AngularPositionModel "${SRCROOT}${OBSERVATIONMODELSDIR}")
target_link_libraries(test_AngularPositionModel tudat_simulation_setup tudat_observation_models tudat_relativity tudat_propagators tudat_acceleration_partials tudat_estimatable_parameters tudat_orbit_determination tudat_aerodynamics tudat_geometric_shapes tudat_gravitation tudat_mission_segments tudat_electro_magnetism tudat_propulsion tudat_ephemerides tudat_numerical_integrators tudat_reference_frames tudat_basic_astrodynamics tudat_input_output tudat_basic_mathematics tudat_propagators tudat_spice_interface cspice ${Boost_LIBRARIES})

    add_executable(test_OneWayRangeModel "${SRCROOT}${OBSERVATIONMODELSDIR}/UnitTests/unitTestOneWayRangeObservationModel.cpp")
    setup_custom_test_program(test_OneWayRangeModel "${SRCROOT}${OBSERVATIONMODELSDIR}")
target_link_libraries(test_OneWayRangeModel tudat_simulation_setup tudat_observation_models tudat_relativity tudat_propagators tudat_acceleration_partials tudat_estimatable_parameters tudat_orbit_determination tudat_aerodynamics tudat_geometric_shapes tudat_gravitation tudat_mission_segments tudat_electro_magnetism tudat_propulsion tudat_ephemerides tudat_numerical_integrators tudat_reference_frames tudat_basic_astrodynamics tudat_input_output tudat_basic_mathematics tudat_propagators tudat_spice_interface cspice ${Boost_LIBRARIES})

    add_executable(test_PositionObservationModel "${SRCROOT}${OBSERVATIONMODELSDIR}/UnitTests/unitTestPositionObservationModel.cpp")
    setup_custom_test_program(test_PositionObservationModel "${SRCROOT}${OBSERVATIONMODELSDIR}")
target_link_libraries(test_PositionObservationModel tudat_simulation_setup tudat_observation_models tudat_relativity tudat_propagators tudat_acceleration_partials tudat_estimatable_parameters tudat_orbit_determination tudat_aerodynamics tudat_geometric_shapes tudat_gravitation tudat_mission_segments tudat_electro_magnetism tudat_propulsion tudat_ephemerides tudat_numerical_integrators tudat_reference_frames tudat_basic_astrodynamics tudat_input_output tudat_basic_mathematics tudat_propagators tudat_spice_interface cspice ${Boost_LIBRARIES})

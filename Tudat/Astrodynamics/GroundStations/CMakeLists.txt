# \file CMakeLists.txt
 #    This is the CMake file for the MissionSegments sub-directory.
 #
 #    Path              : /Astrodynamics/MissionSegments/
 #    Version           : 3
 #    Check status      : Checked
 #
 #    Author            : S.M. Persson
 #    Affiliation       : McGill University, Department of Mechanical Engineering
 #    E-mail address    : mikael.s.persson@gmail.com
 #
 #    Author            : K. Kumar
 #    Affiliation       : Delft University of Technology
 #    E-mail address    : K.Kumar@tudelft.nl
 #
 #    Author            : S.M. Persson
 #    Affiliation       : McGill University, Department of Mechanical Engineering
 #    E-mail address    : mikael.s.persson@gmail.com
 #
 #    Date created      : 20 August, 2011
 #    Last modified     : 26 October, 2011
 #
 #    References
 #
 #    Notes
 #
 #    Copyright (c) 2010-2011 Delft University of Technology.
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

# Set the source files.
set(GROUND_STATIONS_SOURCES
  "${SRCROOT}${GROUNDSTATIONSDIR}/groundStation.cpp"
  "${SRCROOT}${GROUNDSTATIONSDIR}/nominalGroundStationState.cpp"
  "${SRCROOT}${GROUNDSTATIONSDIR}/pointingAnglesCalculator.cpp"
  "${SRCROOT}${GROUNDSTATIONSDIR}/basicTidalBodyDeformation.cpp"
  "${SRCROOT}${GROUNDSTATIONSDIR}/bodyDeformationModel.cpp"
)

# Set the header files.
set(GROUND_STATIONS_HEADERS
  "${SRCROOT}${GROUNDSTATIONSDIR}/groundStation.h"
  "${SRCROOT}${GROUNDSTATIONSDIR}/nominalGroundStationState.h"
  "${SRCROOT}${GROUNDSTATIONSDIR}/pointingAnglesCalculator.h"
  "${SRCROOT}${GROUNDSTATIONSDIR}/basicTidalBodyDeformation.h"
  "${SRCROOT}${GROUNDSTATIONSDIR}/bodyDeformationModel.h"
)

# Add static libraries.
add_library(tudat_ground_stations STATIC ${GROUND_STATIONS_SOURCES} ${GROUND_STATIONS_HEADERS})
setup_tudat_library_target(tudat_ground_stations "${SRCROOT}{GROUNDSTATIONSDIR}")

#add_executable(test_GroundStations "${SRCROOT}${GROUNDSTATIONSDIR}/UnitTests/unitTestGroundStations.cpp")
#setup_custom_test_program(test_GroundStations "${SRCROOT}${GROUNDSTATIONSDIR}")
#target_link_libraries(test_GroundStations tudat_simulation_setup tudat_propagators tudat_aerodynamics tudat_geometric_shapes tudat_gravitation tudat_mission_segments tudat_electro_magnetism tudat_ephemerides tudat_numerical_integrators tudat_reference_frames tudat_basic_astrodynamics tudat_input_output tudat_basic_mathematics tudat_propagators tudat_spice_interface cspice ${Boost_LIBRARIES})

add_executable(test_PointingAnglesCalculator "${SRCROOT}${GROUNDSTATIONSDIR}/UnitTests/unitTestPointingAnglesCalculator.cpp")
setup_custom_test_program(test_PointingAnglesCalculator "${SRCROOT}${GROUNDSTATIONSDIR}")
target_link_libraries(test_PointingAnglesCalculator tudat_simulation_setup tudat_ground_stations tudat_propagators tudat_aerodynamics tudat_geometric_shapes tudat_gravitation tudat_mission_segments tudat_electro_magnetism tudat_ephemerides tudat_numerical_integrators tudat_reference_frames tudat_basic_astrodynamics tudat_input_output tudat_basic_mathematics tudat_propagators tudat_spice_interface cspice ${Boost_LIBRARIES})


#add_executable(test_GroundStationState "${SRCROOT}${GROUNDSTATIONSDIR}/UnitTests/unitTestGroundStationState.cpp")
#setup_custom_test_program(test_GroundStationState "${SRCROOT}${GROUNDSTATIONSDIR}")
#target_link_libraries(test_GroundStationState tudat_simulation_setup tudat_propagators tudat_aerodynamics tudat_geometric_shapes tudat_gravitation tudat_mission_segments tudat_electro_magnetism tudat_ephemerides tudat_numerical_integrators tudat_reference_frames tudat_basic_astrodynamics tudat_input_output tudat_basic_mathematics tudat_propagators tudat_spice_interface cspice ${Boost_LIBRARIES})


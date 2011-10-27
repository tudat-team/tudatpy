# \file CMakeLists.txt
 #    This is the CMake file for the Astrodynamics sub-directory.
 #
 #    Path              : /Astrodynamics/
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

# Define the main sub-directories.
set(BODIESDIR "${ASTRODYNAMICSDIR}/Bodies")
set(ENVIRONMENTMODELSDIR "${ASTRODYNAMICSDIR}/EnvironmentModels")
set(FORCEMODELSDIR "${ASTRODYNAMICSDIR}/ForceModels")
set(MISSIONSEGMENTSDIR "${ASTRODYNAMICSDIR}/MissionSegments")
set(MOMENTMODELSDIR "${ASTRODYNAMICSDIR}/MomentModels")
set(PROPAGATORSDIR "${ASTRODYNAMICSDIR}/Propagators")
set(REFERENCEFRAMESDIR "${ASTRODYNAMICSDIR}/ReferenceFrames")
set(STATESDIR "${ASTRODYNAMICSDIR}/States")

# Set the source files.
set(ASTRODYNAMICS_SOURCES 
  "${SRCROOT}${ASTRODYNAMICSDIR}/physicalConstants.cpp"
)

# Set the header files.
set(ASTRODYNAMICS_HEADERS 
  "${SRCROOT}${ASTRODYNAMICSDIR}/physicalConstants.h"
)

# Add subdirectories.
add_subdirectory("${SRCROOT}${BODIESDIR}")
add_subdirectory("${SRCROOT}${ENVIRONMENTMODELSDIR}")
add_subdirectory("${SRCROOT}${FORCEMODELSDIR}")
add_subdirectory("${SRCROOT}${MISSIONSEGMENTSDIR}")
add_subdirectory("${SRCROOT}${MOMENTMODELSDIR}")
add_subdirectory("${SRCROOT}${PROPAGATORSDIR}")
add_subdirectory("${SRCROOT}${REFERENCEFRAMESDIR}")
add_subdirectory("${SRCROOT}${STATESDIR}")

# Get target properties for static libraries.
get_target_property(BODIESSOURCES tudat_bodies SOURCES)
get_target_property(ENVIRONMENTMODELSSOURCES tudat_environment_models SOURCES)
get_target_property(FORCEMODELSSOURCES tudat_force_models SOURCES)
get_target_property(MISSIONSEGMENTSSOURCES tudat_mission_segments SOURCES)
get_target_property(MOMENTMODELSSOURCES tudat_moment_models SOURCES)
get_target_property(PROPAGATORSSOURCES tudat_propagators SOURCES)
get_target_property(REFERENCEFRAMESSOURCES tudat_reference_frames SOURCES)
get_target_property(STATESSOURCES tudat_states SOURCES)

# Add static libraries.
add_library(tudat_astrodynamics STATIC ${ASTRODYNAMICS_SOURCES} ${ASTRODYNAMICS_HEADERS} ${BODIESSOURCES} ${ENVIRONMENTMODELSSOURCES} ${FORCEMODELSSOURCES} ${MISSIONSEGMENTSSOURCES} ${MOMENTMODELSSOURCE} ${PROPAGATORSSOURCES} ${REFERENCEFRAMESSOURCES} ${STATESSOURCES})
setup_tudat_library_target(tudat_astrodynamics "${SRCROOT}${ASTRODYNAMICSDIR}")
setup_headers("${ASTRODYNAMICS_HEADERS}" "${ASTRODYNAMICSDIR}")

# Add unit tests.
add_executable(test_PhysicalConstants "${SRCROOT}${ASTRODYNAMICSDIR}/unitTestPhysicalConstants.cpp")
setup_custom_test_program(test_PhysicalConstants "${SRCROOT}${ASTRODYNAMICSDIR}")
target_link_libraries(test_PhysicalConstants tudat_input tudat_output tudat_mathematics tudat_astrodynamics tudat_basics)

# End of file.

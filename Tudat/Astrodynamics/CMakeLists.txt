 #    Copyright (c) 2010-2015, Delft University of Technology
 #    All rights reserved.
 #
 #    Redistribution and use in source and binary forms, with or without modification, are
 #    permitted provided that the following conditions are met:
 #      - Redistributions of source code must retain the above copyright notice, this list of
 #        conditions and the following disclaimer.
 #      - Redistributions in binary form must reproduce the above copyright notice, this list of
 #        conditions and the following disclaimer in the documentation and/or other materials
 #        provided with the distribution.
 #      - Neither the name of the Delft University of Technology nor the names of its contributors
 #        may be used to endorse or promote products derived from this software without specific
 #        prior written permission.
 #
 #    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 #    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 #    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 #    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 #    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 #    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 #    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 #    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 #    OF THE POSSIBILITY OF SUCH DAMAGE.
 #
 #    Changelog
 #      YYMMDD    Author            Comment
 #      110820    S.M. Persson      File created.
 #      111025    K. Kumar          Adapted file to work with Revision 194.
 #      111026    K. Kumar          Adapted file so all headers show in project tree in Qt Creator.
 #      120202    K. Kumar          Removed astrodynamicsFunctions.h/.cpp and unit tests.
 #      121123    D. Dirkx          Added ElectroMagnetism sub-directory.
 #
 #    References
 #
 #    Notes
 #

# Define the main sub-directories.
set(AERODYNAMICSDIR "${ASTRODYNAMICSDIR}/Aerodynamics")
set(BASICASTRODYNAMICSDIR "${ASTRODYNAMICSDIR}/BasicAstrodynamics")
set(ELECTROMAGNETISMDIR "${ASTRODYNAMICSDIR}/ElectroMagnetism")
set(EPHEMERIDESDIR "${ASTRODYNAMICSDIR}/Ephemerides")
set(GRAVITATIONDIR "${ASTRODYNAMICSDIR}/Gravitation")
set(MISSIONSEGMENTSDIR "${ASTRODYNAMICSDIR}/MissionSegments")
set(REFERENCEFRAMESDIR "${ASTRODYNAMICSDIR}/ReferenceFrames")
set(OBSERVATIONMODELSDIR "${ASTRODYNAMICSDIR}/ObservationModels")
set(STATEDERIVATIVEMODELSDIR "${ASTRODYNAMICSDIR}/StateDerivativeModels")
set(PROPAGATORSDIR "${ASTRODYNAMICSDIR}/Propagators")

# Add subdirectories.
add_subdirectory("${SRCROOT}${AERODYNAMICSDIR}")
add_subdirectory("${SRCROOT}${BASICASTRODYNAMICSDIR}")
add_subdirectory("${SRCROOT}${ELECTROMAGNETISMDIR}")
add_subdirectory("${SRCROOT}${EPHEMERIDESDIR}")
add_subdirectory("${SRCROOT}${GRAVITATIONDIR}")
add_subdirectory("${SRCROOT}${MISSIONSEGMENTSDIR}")
add_subdirectory("${SRCROOT}${REFERENCEFRAMESDIR}")
add_subdirectory("${SRCROOT}${OBSERVATIONMODELSDIR}")
add_subdirectory("${SRCROOT}${STATEDERIVATIVEMODELSDIR}")
add_subdirectory("${SRCROOT}${PROPAGATORSDIR}")

# Get target properties for static libraries.
get_target_property(AERODYNAMICSSOURCES tudat_aerodynamics SOURCES)
get_target_property(BASICASTRODYNAMICSSOURCES tudat_basic_astrodynamics SOURCES)
get_target_property(ELECTROMAGNETISMSOURCES tudat_electro_magnetism SOURCES)
get_target_property(EPHEMERIDESSOURCES tudat_ephemerides SOURCES)
get_target_property(GRAVITATIONSOURCES tudat_gravitation SOURCES)
get_target_property(MISSIONSEGMENTSSOURCES tudat_mission_segments SOURCES)
get_target_property(REFERENCEFRAMESSOURCES tudat_reference_frames SOURCES)
get_target_property(OBSERVATIONMODELSDIR tudat_observation_models SOURCES)
get_target_property(STATEDERIVATIVEMODELSSOURCES tudat_state_derivative_models SOURCES)
get_target_property(PROPAGATORSDIR tudat_propagators SOURCES)

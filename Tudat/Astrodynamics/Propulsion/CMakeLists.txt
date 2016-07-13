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
 #	120203    K. Kumar	    Adapted for new BasicAstrodynamics package.
 #      120216    K. Kumar          Removed propagator files.
 #      120822    P. Musegaas       Added unit test for mean to hyperbolic eccentric anomaly
 #                                  conversion.
 #      130218    D. Dirkx          Added unit test for Julian date conversions.
 #      130214    R.C.A. Boon       Added modified equinoctial element files (.h/.cpp, unit test).
 #      130301    D. Dirkx          Added unit test for geodetic coordinate conversions.
 #      160406    M. Van den Broeck Added Unified State Model element files (.h/.cpp, unit test).
 #
 #    References
 #
 #    Notes
 #

# Set the source files.
set(PROPULSION_SOURCES
  "${SRCROOT}${PROPULSIONDIR}/thrustGuidance.cpp"
  "${SRCROOT}${PROPULSIONDIR}/thrustFunctions.cpp"
)

# Set the header files.
set(PROPULSION_HEADERS
  "${SRCROOT}${PROPULSIONDIR}/thrustAccelerationModel.h"
  "${SRCROOT}${PROPULSIONDIR}/thrustGuidance.h"
  "${SRCROOT}${PROPULSIONDIR}/thrustFunctions.h"
  "${SRCROOT}${PROPULSIONDIR}/massRateFromThrust.h"
)

# Add static libraries.
add_library(tudat_propulsion STATIC ${PROPULSION_SOURCES} ${PROPULSION_HEADERS})
setup_tudat_library_target(tudat_propulsion "${SRCROOT}${PROPULSIONDIR}")

# Add unit tests.
set(PROPULSION_UNITTESTS
  "${SRCROOT}${PROPULSIONDIR}/UnitTests/unitTestBasicAstrodynamics.cpp"
  "${SRCROOT}${PROPULSIONDIR}/UnitTests/unitTestAstrodynamicsFunctions.cpp"
  "${SRCROOT}${PROPULSIONDIR}/UnitTests/unitTestOrbitalElementConversions.cpp"
  "${SRCROOT}${PROPULSIONDIR}/UnitTests/unitTestPhysicalConstants.cpp"
  "${SRCROOT}${PROPULSIONDIR}/UnitTests/unitTestUnitConversions.cpp"
  "${SRCROOT}${PROPULSIONDIR}/UnitTests/unitTestUnifiedStateModelElementConversions.cpp"
)

add_executable(test_ThrustAcceleration "${SRCROOT}${PROPULSIONDIR}/UnitTests/unitTestThrustAcceleration.cpp")
setup_custom_test_program(test_ThrustAcceleration "${SRCROOT}${PROPULSIONDIR}")
target_link_libraries(test_ThrustAcceleration tudat_simulation_setup tudat_propagators tudat_acceleration_partials tudat_estimatable_parameters tudat_orbit_determination tudat_aerodynamics tudat_geometric_shapes tudat_gravitation tudat_mission_segments tudat_electro_magnetism tudat_propulsion tudat_ephemerides tudat_numerical_integrators tudat_reference_frames tudat_basic_astrodynamics tudat_input_output tudat_basic_mathematics tudat_propagators tudat_spice_interface cspice ${Boost_LIBRARIES})

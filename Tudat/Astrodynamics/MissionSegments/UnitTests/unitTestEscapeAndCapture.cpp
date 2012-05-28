/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110131    E. Iorfida        First creation of the code.
 *      110202    E. Iorfida        Modified structure of the code, unique
 *                                  base class for launch and capture paths.
 *      110208    E. Iorfida        Deleted execute( ) function. Modified getDeltaV into
 *                                  computeDeltaV.
 *      110214    E. Iorfida        Code updated with the modification made in .h/.cpp files
 *                                  about radius of central body.
 *      110627    K. Kumar          Updated to use new predefined planets code.
 *      120416    T. Secretin       Boostified unit test.
 *
 *    References        :
 *      Mengali, G., Quarta, A.A. Fondamenti di Meccanica del volo Spaziale,
 *          Edizioni PLUS.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/Bodies/planet.h"
#include "Tudat/Astrodynamics/MissionSegments/capturePhase.h"
#include "Tudat/Astrodynamics/MissionSegments/escapeAndCapture.h"
#include "Tudat/Astrodynamics/MissionSegments/escapePhase.h"

namespace tudat
{
namespace unit_tests
{

//! Test patched conics implementation.
BOOST_AUTO_TEST_SUITE( test_escape_and_capture )

//! Test delta-V computation for the escape phase.
BOOST_AUTO_TEST_CASE( testDeltaVEscape )
{
    // Set tolerance.
    const double tolerance = 1.0e-1;

    // Expected test result based on example 11.4 page 334, "Fondamenti di
    // Meccanica del volo Spaziale", G. Mengali, A. A. Quarta.
    const double expectedDeltaVEscape = 3.5244e3;

    // Set test case.
    using astrodynamics::mission_segments::EscapePhase;
    EscapePhase escapePhaseTest;

    // Central body parameters.
    bodies::Planet predefinedEarth;
    predefinedEarth.setPredefinedPlanetSettings( bodies::Planet::earth );
    
    // Set launch conditions.
    escapePhaseTest.setCentralGravityField( predefinedEarth.getGravityFieldModel( ) );
    escapePhaseTest.setParkingOrbitRadius( 6371.0e3 );
    escapePhaseTest.setPeriapsisAltitude( 629.0e3 );
    escapePhaseTest.setEccentricity( 0.0 );
    escapePhaseTest.setHyperbolicExcessSpeed( 2.9444e3 );

    // Compute delta-V of escape phase.
    const double deltaVEscape_ = escapePhaseTest.computeDeltaV( );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaVEscape, deltaVEscape_, tolerance );
}

//! Test delta-V computation for the capture phase.
BOOST_AUTO_TEST_CASE( testDeltaVCapture )
{
    // Set tolerance.
    const double tolerance = 1.0e-1;

    // Expected test result based on example 11.4 page 334, "Fondamenti di
    // Meccanica del volo Spaziale", G. Mengali, A. A. Quarta.
    const double expectedDeltaVCapture = 1.9425e3;

    // Set test case.
    using astrodynamics::mission_segments::CapturePhase;
    CapturePhase capturePhaseTest;

    // Central body parameters.
    bodies::Planet predefinedMars;
    predefinedMars.setPredefinedPlanetSettings( bodies::Planet::mars );

    // Set capture conditions.
    capturePhaseTest.setCentralGravityField( predefinedMars.getGravityFieldModel( ) );
    capturePhaseTest.setParkingOrbitRadius( 3389.0e3 );
    capturePhaseTest.setPeriapsisAltitude( 2611.0e3 );
    capturePhaseTest.setEccentricity( 0.0 );
    capturePhaseTest.setHyperbolicExcessSpeed( 2.6486e3 );

    // Compute delta-V of capture phase.
    const double deltaVCapture_ = capturePhaseTest.computeDeltaV( );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaVCapture, deltaVCapture_, tolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

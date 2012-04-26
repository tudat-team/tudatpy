/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
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

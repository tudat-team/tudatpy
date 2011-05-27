/*! \file unitTestEscapeAndCapture.cpp
 *    Source file of unit test for escapePhase method.
 *    The unit test is based on the example 11.4 page 334, "Fondamenti di
 *    Meccanica del volo Spaziale", G. Mengali, A. A. Quarta.
 *
 *    Path              : /Astrodynamics/MissionSegments/EscapeAndCapture/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author            : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 31 January, 2011
 *    Last modified     : 14 February, 2011
 *
 *    References        :
 *      Mengali, G., Quarta, A.A. Fondamenti di Meccanica del volo Spaziale,
 *          Edizioni PLUS.
 *
 *    Notes
 *      Test runs code and verifies result against expected value.
 *      If the tested code is erroneous, the test function returns a boolean
 *      true; if the code is correct, the function returns a boolean false.
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 *      110208    E. Iorfida        Deleted execute() function. Modified
 *                                  getDeltaV into computeDeltaV.
 *      110214    E. Iorfida        Code updated with the modification made
 *                                  in .h/.cpp files about radius of central
 *                                  body.
 */

// Include statements.
#include "unitTestEscapeAndCapture.h"
#include "sphereSegment.h"

// Using declarations.
using std::cerr;
using std::endl;
using predefined_celestial_bodies::createPredefinedCelestialBody;
using mathematics::computeAbsoluteValue;

//! Namespace for all unit tests.
namespace unit_tests
{

//! Test patched conics implementation.
bool testEscapeAndCapture( )
{
    // Test result initialised to false.
    bool isEscapeAndCaptureErroneous = false;

    // Set tolerance.
    double tolerance = 1.0e-1;

    // Expected test result based on example 11.4 page 334, "Fondamenti di
    // Meccanica del volo Spaziale", G. Mengali, A. A. Quarta.
    double expectedDeltaVEscape = 3.5244e3;
    double expectedDeltaVCapture = 1.9425e3;

    // Set test case.
    EscapePhase* pointerToEscapePhaseTest = new EscapePhase;
    CapturePhase* pointerToCapturePhaseTest = new CapturePhase;

    // Central bodies parameters.
    // Central body at launch phase.
    CelestialBody* pointerToEarth = new CelestialBody;
    pointerToEarth = createPredefinedCelestialBody(
            predefined_celestial_bodies::earth );
    SphereSegment* pointerToEarthSphere = new SphereSegment;
    pointerToEarth->setShapeModel( pointerToEarthSphere );

    // Central body at capture phase.
    CelestialBody* pointerToMars = new CelestialBody;
    pointerToMars = createPredefinedCelestialBody(
            predefined_celestial_bodies::mars );
    SphereSegment* pointerToMarsSphere = new SphereSegment;
    pointerToMars->setShapeModel( pointerToMarsSphere );

    // Set launch conditions.
    pointerToEscapePhaseTest->setCentralBody( pointerToEarth );
    pointerToEarthSphere->setRadius( 6371.0e3 );
    pointerToEscapePhaseTest->setPeriapsisAltitude( 629.0e3 );
    pointerToEscapePhaseTest->setEccentricity( 0.0 );
    pointerToEscapePhaseTest->setHyperbolicExcessSpeed( 2.9444e3 );

    // Set capture conditions.
    pointerToCapturePhaseTest->setCentralBody( pointerToMars );
    pointerToMarsSphere->setRadius( 3389.0e3 );
    pointerToCapturePhaseTest->setPeriapsisAltitude( 2611.0e3 );
    pointerToCapturePhaseTest->setEccentricity( 0.0 );
    pointerToCapturePhaseTest->setHyperbolicExcessSpeed( 2.6486e3 );

    // Execute patched conic implementation.
    pointerToEscapePhaseTest->computeDeltaV( );
    pointerToCapturePhaseTest->computeDeltaV( );

    // Define delta-V of escape/capture phase.
    double deltaVEscape_ = pointerToEscapePhaseTest->computeDeltaV( );
    double deltaVCapture_ = pointerToCapturePhaseTest->computeDeltaV( );

    // Set test result to true if the test does not match the expected result.
    if ( mathematics::computeAbsoluteValue(
         deltaVEscape_ - expectedDeltaVEscape ) >= tolerance ||
         mathematics::computeAbsoluteValue(
         deltaVCapture_ - expectedDeltaVCapture ) >= tolerance )
    {
        isEscapeAndCaptureErroneous = true;

        if ( mathematics::computeAbsoluteValue(
             deltaVEscape_ - expectedDeltaVEscape ) >= tolerance )
        {
            cerr << "The computed value of the delta-V of the launch phase ("
                 << deltaVEscape_
                 << ") does not match "
                 << "the expected solution (" << expectedDeltaVEscape << ")."
                 << endl;
            cerr << "The difference is: "
                 << mathematics::computeAbsoluteValue(
                    deltaVEscape_ - expectedDeltaVEscape )
                 << endl;
        }
        else if ( mathematics::computeAbsoluteValue(
                  deltaVCapture_ - expectedDeltaVCapture ) >= tolerance )
        {
            cerr << "The computed value of the delta-V of the capture phase ("
                 << deltaVCapture_
                 << ") does not match "
                 << "the expected solution (" << expectedDeltaVCapture << ")."
                 << endl;
            cerr << "The difference is: "
                 << mathematics::computeAbsoluteValue(
                    deltaVCapture_ - expectedDeltaVCapture )
                 << endl;
        }
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isEscapeAndCaptureErroneous;
}

}

// End of file.

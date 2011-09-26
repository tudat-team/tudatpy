/*! \file unitTestReferenceFrameTransformations.cpp
 *    This file contains the implementation of the reference frame transformation
 *    unit test included in Tudat.
 *
 *    Path              : /Astrodynamics/ReferenceFrames/
 *    Version           : 6
 *    Check status      : Checked
 *
 *    Checker           : F. M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 19 May, 2011
 *    Last modified     : 08 August, 2011
 *
 *    References
 *
 *    Notes
 *      The reference frame definitions/abbreviations can be found in the file
 *      referenceFrameTransformations.h.
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      110519    F.M. Engelen      First creation of code.
 *      110628    K. Kumar          Minor comment and layout modifications.
 *      110701    K. Kumar          Updated unit tests to check relative error;
 *                                  Updated file path.
 *      110718    F.M. Engelen      Took out falacy in test (only checking the norm) and
 *                                  added ItoE and EtoI transformation.
 *      110726    K. Kumar          Minor modifications; updated relative error
 *                                  wrt to norm.
 *      110808    F.M. Engelen      Updated with better tests, changed test for vertical frame.
 */

// Include statements.
#include <cmath>
#include "unitTestReferenceFrameTransformations.h"

//using directives
using std::cerr;
using std::endl;
using std::fabs;
using mathematics::MACHINE_PRECISION_DOUBLES;
using unit_conversions::convertDegreesToRadians;

//! Namespace for all unit tests.
namespace unit_tests
{

//! Test reference frame transformations.
bool testReferenceFrameTransformations( )
{
    // Summary of tests.
    // Test 1: Test inertial to rotating planetocentric frame transformation.
    // Test 2: Check if the transformed matrix of Test 1 is also correct.
    // Test 3: Same test as Test 1 for the transformation quaternion.
    // Test 4: Same test as Test 2 for the transformation quaternion.
    // Test 5: Test airspeed-based aerodynamic to body frame transformation.
    // Test 6: Test planetocentric to local vertical frame transformation quaternion.
    // Test 7: Check if the transformed matrix of Test 6 is also correct.

    // Test 1: Test Inertial to rotating planetocentric frame transformation.

    // Initialize the test boolean.
    bool isFrameTransformationErroneous = false;

    // Initialize initial location vector in inertial frame.
    Vector3d startLocation;
    startLocation( 0 ) = 4.0;
    startLocation( 1 ) = 1.0;
    startLocation( 2 ) = 5.0;

    double horizontalStartLocationSize = sqrt( pow( startLocation( 0 ), 2.0 ) + pow( startLocation( 1 ), 2.0 ) );

    // Declare and initialize angle between vector and XR-axis.
    double startAngle = atan2( startLocation( 1 ) , startLocation( 0 ) );

    // Rotate by 10 degrees around the positive Z-axis
    double angleInTime = convertDegreesToRadians( 10.0 );

    // Declare and initialize the angle between the XI- and XR-axis.
    double endAngle = startAngle - angleInTime;

    // Declare the expected location of the point in the planetocentric reference frame.
    Vector3d expectedLocation;
    expectedLocation( 0 ) = horizontalStartLocationSize * cos( endAngle );
    expectedLocation( 1 ) = horizontalStartLocationSize * sin( endAngle );
    expectedLocation( 2 ) = startLocation( 2 );

    // Compute location of the point in the totating frame subject to the transformation matrix.
    Vector3d transformedLocation = reference_frame_transformations::
            getInertialToPlanetocentricFrameTransformationMatrix( angleInTime ) * startLocation;

    // Compute the error in the calculation.
    Vector3d absoluteNumericalError =  transformedLocation - expectedLocation;
    double relativeNumericalError;

    // Check if relative errors are too large and output cerr statements if necessary.
    for ( int i = 0; i < absoluteNumericalError.size( ); i++ )
    {
        // Compute relative error of given component.
        relativeNumericalError = absoluteNumericalError( i ) / expectedLocation.norm( );

        if ( fabs( relativeNumericalError ) > MACHINE_PRECISION_DOUBLES )
        {
            isFrameTransformationErroneous = true;

            cerr << "Frame transformation Test 1." << i << " failed " << endl;
            cerr << transformedLocation( i ) << endl;
            cerr << expectedLocation( i ) << endl;
        }
    }

    // Test 2: Check if the transformed matrix of Test 1 is also correct.
    // Compute the error in the calculation.
    absoluteNumericalError =
            reference_frame_transformations::
            getRotatingPlanetocentricToInertialFrameTransformationMatrix( angleInTime ) *
            reference_frame_transformations::
            getInertialToPlanetocentricFrameTransformationMatrix( angleInTime )
            * startLocation
            - startLocation;

    // Check if relative errors are too large and output cerr statements if necessary.
    for ( int i = 0; i < absoluteNumericalError.size( ); i++ )
    {
        // Compute relative error of given component.
        relativeNumericalError = absoluteNumericalError( i ) / expectedLocation.norm( );

        if ( fabs( relativeNumericalError ) > MACHINE_PRECISION_DOUBLES )
        {
            isFrameTransformationErroneous = true;

            cerr << "Frame transformation Test 2." << i << " failed " << endl;
            cerr << transformedLocation( i ) << endl;
            cerr << expectedLocation( i ) << endl;
        }
    }

    // Test 3: Same test as Test 1 for the transformation quaternion.
    // Compute location of the point in the Rotating frame subject to the transformation matrix.
    transformedLocation = reference_frame_transformations::
            getInertialToPlanetocentricFrameTransformationQuaternion( angleInTime ) * startLocation;

    // Compute the error in the calculation.
    absoluteNumericalError = transformedLocation - expectedLocation;

    // Check if relative errors are too large and output cerr statements if necessary.
    for ( int i = 0; i < absoluteNumericalError.size( ); i++ )
    {
        // Compute relative error of given component.
        relativeNumericalError = absoluteNumericalError( i ) / expectedLocation.norm( );

        if ( fabs( relativeNumericalError ) > MACHINE_PRECISION_DOUBLES )
        {
            isFrameTransformationErroneous = true;

            cerr << "Frame transformation Test 3." << i << " failed " << endl;
            cerr << transformedLocation( i ) << endl;
            cerr << expectedLocation( i ) << endl;
        }
    }

    // Test 4: Same test as Test 2 for the transformation quaternion.
    // Compute the error in the calculation.
    absoluteNumericalError =
            reference_frame_transformations::
            getRotatingPlanetocentricToInertialFrameTransformationQuaternion( angleInTime ) *
            reference_frame_transformations::
            getInertialToPlanetocentricFrameTransformationQuaternion( angleInTime )
            * startLocation
            - startLocation;

    // Check if relative errors are too large and output cerr statements if necessary.
    for ( int i = 0; i < absoluteNumericalError.size( ); i++ )
    {
        // Compute relative error of given component.
        relativeNumericalError = absoluteNumericalError( i ) / expectedLocation.norm( );

        if ( fabs( relativeNumericalError ) > MACHINE_PRECISION_DOUBLES )
        {
            isFrameTransformationErroneous = true;

            cerr << "Frame transformation Test 4." << i << " failed " << endl;
            cerr << transformedLocation( i ) << endl;
            cerr << expectedLocation( i ) << endl;
        }
    }

    // Test 5: Test airspeed-Based Aerodynamic to body frame transformation.
    // Declare and initialize the start location and angles.
    startLocation( 0 ) = -10.0;
    startLocation( 1 ) = 0.0;
    startLocation( 2 ) = 0.0;
    double angleOfAttack = convertDegreesToRadians( 45.0 );
    double angleOfSideslip = convertDegreesToRadians( 60.0 );;

    // Compute expected location.
    // As there is only an angle of attack, the following simplified equations can be used.
    expectedLocation( 0 ) = startLocation( 0 ) * cos( angleOfSideslip ) * cos( angleOfAttack );
    expectedLocation( 1 ) = startLocation( 0 ) * sin( angleOfSideslip );
    expectedLocation( 2 ) = startLocation( 0 ) * cos( angleOfSideslip ) * sin( angleOfAttack );

    // Transform the vector from the AA to the B frame.
    transformedLocation = reference_frame_transformations::
            getAirspeedBasedAerodynamicToBodyFrameTransformationMatrix( angleOfAttack,
                                                                        angleOfSideslip )
            * startLocation;

    // Calculate error.
    absoluteNumericalError = transformedLocation - expectedLocation;

    // Check if relative errors are too large and output cerr statements if necessary.
    for ( int i = 0; i < absoluteNumericalError.size( ); i++ )
    {
        // Compute relative error of given component.
        relativeNumericalError = absoluteNumericalError( i ) / expectedLocation.norm( );

        if ( fabs( relativeNumericalError ) > 1.0e-14 )
        {
            isFrameTransformationErroneous = true;

            cerr << "Frame transformation Test 5." << i << " failed " << endl;
            cerr << transformedLocation( i ) << endl;
            cerr << expectedLocation( i ) << endl;
            cerr << relativeNumericalError << " " << MACHINE_PRECISION_DOUBLES << endl;
        }
    }


    // Test 6: Test Rotating planetocentric to local vertical frame transformation quaternion.

    // Initialize initial location vector.
    startLocation( 0 ) = 10.0;
    startLocation( 1 ) = 5.0;
    startLocation( 2 ) = 2.0;

    // Declare rotation angle for planet and set to 90 degrees.
    double longitude = convertDegreesToRadians( 60.0 );

    // Declare latitude and set to 45 degrees.
    double latitude = convertDegreesToRadians( 20.0 );

    // Declare the expected location of the point in the planet reference frame.
    expectedLocation( 0 ) = -cos( longitude ) * sin( latitude ) * startLocation( 0 ) -
                            sin( longitude ) * sin( latitude ) * startLocation( 1 )  +
                            cos( latitude ) * startLocation( 2 );
    expectedLocation( 1 ) = -sin( longitude ) * startLocation( 0 ) +
                            cos( longitude ) * startLocation( 1 ) +
                            0.0;
    expectedLocation( 2 ) = -cos( longitude ) * cos( latitude ) * startLocation( 0 ) -
                            sin( longitude ) * cos( latitude ) * startLocation( 1 )  -
                            sin( latitude ) * startLocation( 2 );

    // Compute location of the point in the inertial frame subject to the transformation matrix.
    transformedLocation = reference_frame_transformations::
            getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion( longitude, latitude )
            * startLocation;

    // Compute the error in the calculation.
    absoluteNumericalError =  transformedLocation - expectedLocation;

    // Check if relative errors are too large and output cerr statements if necessary.
    for ( unsigned int i = 0; i < 3; i++ )
    {
        // Compute relative error of given component.
        relativeNumericalError = absoluteNumericalError( i );

        if ( fabs( relativeNumericalError ) > 1.0e-14 )
        {
            isFrameTransformationErroneous = true;

            cerr << "Frame transformation Test 6." << i << " failed " << endl;
            cerr << transformedLocation( i ) << endl;
            cerr << expectedLocation( i ) << endl;
            cerr << relativeNumericalError << "  " << MACHINE_PRECISION_DOUBLES << endl;
        }
    }

    // Test 7: Check if the transformed matrix of Test 6 is also correct.
    // Compute the error in the calculation.
    absoluteNumericalError =
            reference_frame_transformations::
            getLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion( longitude, latitude )
            * reference_frame_transformations::
            getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion( longitude, latitude )
            * startLocation - startLocation;

    // Compute relative error in first component.
    relativeNumericalError = absoluteNumericalError( 0 );

    // Check if relative error is too large and output cerr statements if necessary.
    if ( fabs( relativeNumericalError ) > 1.0e-14 )
    {
        isFrameTransformationErroneous = true;

        cerr << "Frame transformation Test 7." << 0 << " failed " << endl;
    }

    // Check if relative error of remaining compoennts is too large and output cerr statements if
    // necessary.
    for ( int i = 1; i< absoluteNumericalError.size( ); i++ )
    {
        relativeNumericalError = absoluteNumericalError( i ) / expectedLocation.norm( );

        if ( fabs( relativeNumericalError ) > 1.0e-14 )
        {
            isFrameTransformationErroneous = true;

            cerr << "Frame transformation Test 7." << i << " failed " << endl;
            cerr << transformedLocation( i ) << endl;
            cerr << expectedLocation( i ) << endl;
            cerr << relativeNumericalError << "  " << MACHINE_PRECISION_DOUBLES << endl;
        }
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isFrameTransformationErroneous;
}

}

// End of file.

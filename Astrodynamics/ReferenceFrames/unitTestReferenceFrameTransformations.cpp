/*! \file referenceFrameTransformations.cpp
 *    This file contains the implementation of the reference frame
 *    transformation unit test included in Tudat.
 *
 *    Path              : /Astrodynamics/ReferenceFrames/
 *    Version           : 3
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
 *    Last modified     : 1 July, 2011
 *
 *    References
 *
 *    Notes
 *    The reference frame definitions/abbreviations can be found in the file
 *      referenceFrameTransformations.h.
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
 *      110519    F.M. Engelen      First creation of code.
 *      110628    K. Kumar          Minor comment and layout modifications.
 *      110701    K. Kumar          Updated unit tests to check relative error;
 *                                  Updated file path.
 */

// Include statements.
#include "unitTestReferenceFrameTransformations.h"

//using directives
using std::cerr;
using std::endl;
using mathematics::computeAbsoluteValue;
using mathematics::MACHINE_PRECISION_DOUBLES;
using unit_conversions::convertDegreesToRadians;

//! Namespace for all unit tests.
namespace unit_tests
{

//! Test reference frame transformations.
bool testReferenceFrameTransformations( )
{
    // Summary of tests.
    // Test 1: Test rotating planetocentric to Inertial frame transformation.
    // Test 2: Check if the transformed matrix of Test 1 is also correct.
    // Test 3: Test airspeed-Based Aerodynamic to body frame transformation.

    // Test 1: Test rotating planetocentric to Inertial frame transformation.

    // Initialize the test boolean.
    bool isFrameTransformationErroneous = false;

    // Initialize initial location vector.
    Vector3d startLocation;
    startLocation( 0 ) = 4.0;
    startLocation( 1 ) = 3.0;
    startLocation( 2 ) = 5.0;

    // Declare and initialize angle between vector and XR-axis.
    double startAngle = atan( 3.0 / 4.0 );

    // Rotate by 10 degrees around the positive Z-axis
    double angleInTime = convertDegreesToRadians( 10.0 );

    // Declare and initialize the angle between the XI- and XR-axis.
    double endAngle = startAngle + angleInTime;

    // Declare the expected location of the point in the Inertial reference frame.
    Vector3d expectedLocation;
    expectedLocation( 0 ) = 5.0 * cos( endAngle );
    expectedLocation( 1 ) = 5.0 * sin( endAngle );
    expectedLocation( 2 ) = 5.0;

    // Compute location of the point in the Inertial frame subject to the transformation matrix.
    Vector3d transformedLocation = reference_frame_transformations::
            getRotatingPlanetocentricToInertialFrameTransformation( angleInTime ) * startLocation;

    // Compute the error in the calculation.
    double absoluteErrorInTransformation = computeAbsoluteValue(
                ( transformedLocation.norm( ) - expectedLocation.norm( ) ) );

    double relativeNumericalError = absoluteErrorInTransformation /
                                    computeAbsoluteValue( transformedLocation.norm( ) );

    // Check if error is too large (factor 3 because the vector is 3 long).
    if ( relativeNumericalError > MACHINE_PRECISION_DOUBLES )
    {
        // Output error statement.
        cerr << "The Rotating planetocentric (R) to Inertial (I) transformation is incorrect."
             << endl;

        isFrameTransformationErroneous = true;
    }

    // Test 2: Check if the transformed matrix of Test 1 is also correct.

    // Compute the error in the calculation.
    absoluteErrorInTransformation = computeAbsoluteValue(
                ( reference_frame_transformations::
                  getRotatingPlanetocentricToInertialFrameTransformation( angleInTime ) *
                  reference_frame_transformations::
                  getInertialToPlanetocentricFrameTransformation( angleInTime )
                  * startLocation ).norm( )
                - startLocation.norm( ) );

    relativeNumericalError = absoluteErrorInTransformation /
                                        computeAbsoluteValue( startLocation.norm( ) );


    // Check if error is too large (factor 3 because the vector is 3 long).
    if ( relativeNumericalError > MACHINE_PRECISION_DOUBLES )
    {
        // Output error statement.
        cerr << "The Inertial (I) to Rotating planetocentric (R) transformation is incorrect."
             << endl;

        isFrameTransformationErroneous = true;
    }

    // Test 3: Test airspeed-Based Aerodynamic to body frame transformation.

    // Declare and initialize the start location and angles.
    startLocation( 0 ) = -10.0;
    startLocation( 1 ) = 0.0;
    startLocation( 2 ) = 0.0;
    double angleOfAttack = convertDegreesToRadians( 10.0 );
    double angleOfSideslip = 0.0 * M_PI;

    // Compute expected location.
    // As there is only an angle of attack, the following simplified equations can be used.
    expectedLocation( 0 ) = startLocation( 0 ) * cos( angleOfAttack );
    expectedLocation( 1 ) = 0.0;
    expectedLocation( 2 ) = startLocation( 0 ) * sin( angleOfAttack );

    // Transform the vector from the AA to the B frame.
    transformedLocation = reference_frame_transformations::
            getAirspeedBasedAerodynamicToBodyFrameTransformation( angleOfAttack, angleOfSideslip )
            * startLocation;

    // Calculate error.
    absoluteErrorInTransformation = computeAbsoluteValue(
                ( transformedLocation.norm( ) - expectedLocation.norm( ) ) );

    relativeNumericalError = absoluteErrorInTransformation /
                                        computeAbsoluteValue( transformedLocation.norm( ) );

    // Check if error is too large (factor 3 because the vector is 3 long).
    if ( relativeNumericalError > ( MACHINE_PRECISION_DOUBLES ) )
    {
        // Output error statement.
        cerr << "The Aerodynamic (airspeed based) (AA) to body reference frame (B) " << endl;
        cerr << "tranformation is incorrect." << endl;

        isFrameTransformationErroneous = true;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isFrameTransformationErroneous;
}

}

// End of file.

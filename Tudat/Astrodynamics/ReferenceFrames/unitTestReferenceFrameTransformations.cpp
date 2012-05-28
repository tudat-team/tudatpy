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
 *      110519    F.M. Engelen      Creation of code.
 *      110628    K. Kumar          Minor comment and layout modifications.
 *      110701    K. Kumar          Updated unit tests to check relative error;
 *                                  Updated file path.
 *      110718    F.M. Engelen      Took out falacy in test (only checking the norm) and
 *                                  added ItoE and EtoI transformation.
 *      110726    K. Kumar          Minor modifications; updated relative error
 *                                  wrt to norm.
 *      110808    F.M. Engelen      Updated with better tests, changed test for vertical frame.
 *
 *    References
 *
 *    The reference frame definitions/abbreviations can be found in the file
 *    referenceFrameTransformations.h.
 *
 */

#include <cmath>
#include <iostream>
#include <limits>

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

//! Test reference frame transformations.
int main( )
{
    // Using declarations.
    using std::atan2;
    using std::cos;
    using std::sin;
    using std::cerr;
    using std::endl;
    using std::fabs;
    using std::pow;
    using std::sqrt;
    using tudat::unit_conversions::convertDegreesToRadians;
    using namespace tudat;

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
    Eigen::Vector3d startLocation;
    startLocation( 0 ) = 4.0;
    startLocation( 1 ) = 1.0;
    startLocation( 2 ) = 5.0;

    double horizontalStartLocationSize = sqrt( pow( startLocation( 0 ), 2.0 )
                                               + pow( startLocation( 1 ), 2.0 ) );

    // Declare and initialize angle between vector and XR-axis.
    double startAngle = atan2( startLocation( 1 ), startLocation( 0 ) );

    // Rotate by 10 degrees around the positive Z-axis
    double angleInTime = convertDegreesToRadians( 10.0 );

    // Declare and initialize the angle between the XI- and XR-axis.
    double endAngle = startAngle - angleInTime;

    // Declare the expected location of the point in the planetocentric reference frame.
    Eigen::Vector3d expectedLocation;
    expectedLocation( 0 ) = horizontalStartLocationSize * cos( endAngle );
    expectedLocation( 1 ) = horizontalStartLocationSize * sin( endAngle );
    expectedLocation( 2 ) = startLocation( 2 );

    // Compute location of the point in the totating frame subject to the transformation matrix.
    Eigen::Vector3d transformedLocation = reference_frame_transformations::
            getInertialToPlanetocentricFrameTransformationMatrix( angleInTime ) * startLocation;

    // Compute the error in the calculation.
    Eigen::Vector3d absoluteNumericalError = transformedLocation - expectedLocation;
    double relativeNumericalError;

    // Check if relative errors are too large and output cerr statements if necessary.
    for ( int i = 0; i < absoluteNumericalError.size( ); i++ )
    {
        // Compute relative error of given component.
        relativeNumericalError = absoluteNumericalError( i ) / expectedLocation.norm( );

        if ( fabs( relativeNumericalError ) > std::numeric_limits< double >::epsilon( ) )
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
            * startLocation - startLocation;

    // Check if relative errors are too large and output cerr statements if necessary.
    for ( int i = 0; i < absoluteNumericalError.size( ); i++ )
    {
        // Compute relative error of given component.
        relativeNumericalError = absoluteNumericalError( i ) / expectedLocation.norm( );

        if ( fabs( relativeNumericalError ) > std::numeric_limits< double >::epsilon( ) )
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
            getInertialToPlanetocentricFrameTransformationQuaternion( angleInTime )
            * startLocation;

    // Compute the error in the calculation.
    absoluteNumericalError = transformedLocation - expectedLocation;

    // Check if relative errors are too large and output cerr statements if necessary.
    for ( int i = 0; i < absoluteNumericalError.size( ); i++ )
    {
        // Compute relative error of given component.
        relativeNumericalError = absoluteNumericalError( i ) / expectedLocation.norm( );

        if ( fabs( relativeNumericalError ) > std::numeric_limits< double >::epsilon( ) )
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
            * startLocation - startLocation;

    // Check if relative errors are too large and output cerr statements if necessary.
    for ( int i = 0; i < absoluteNumericalError.size( ); i++ )
    {
        // Compute relative error of given component.
        relativeNumericalError = absoluteNumericalError( i ) / expectedLocation.norm( );

        if ( fabs( relativeNumericalError ) > std::numeric_limits< double >::epsilon( ) )
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
            cerr << relativeNumericalError << " " << std::numeric_limits< double >::epsilon( )
                 << endl;
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
            getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion( longitude,
                                                                                   latitude )
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
            cerr << relativeNumericalError << "  " << std::numeric_limits< double >::epsilon( )
                 << endl;
        }
    }

    // Test 7: Check if the transformed matrix of Test 6 is also correct.
    // Compute the error in the calculation.
    absoluteNumericalError =
            reference_frame_transformations::
            getLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion( longitude,
                                                                                   latitude )
            * reference_frame_transformations::
            getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion( longitude,
                                                                                   latitude )
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
            cerr << relativeNumericalError << "  " << std::numeric_limits< double >::epsilon( )
                 << endl;
        }
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    if ( isFrameTransformationErroneous )
    {
        cerr << "referenceFrameTransformations failed!" << endl;
    }

    return isFrameTransformationErroneous;
}

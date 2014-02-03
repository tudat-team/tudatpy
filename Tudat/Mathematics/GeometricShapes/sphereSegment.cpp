/*    Copyright (c) 2010-2014, Delft University of Technology
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
 *      100910    D. Dirkx          First version of file.
 *      100915    D. Dirkx          Modified to correct comments, 80-lines
 *                                  rule, etc.
 *      100928    D. Dirkx          Modifications following first checking
 *                                  iteration.
 *      100929    D. Dirkx          Creation of separate file for class.
 *      101125    D. Dirkx          Update of class, better get and set
 *                                  functions.
 *      110208    K. Kumar          Updated file header; Doxygen comments
 *                                  corrected; minor changes to functions.
 *      110209    D. Dirkx          Minor changes.
 *      110209    K. Kumar          Minor changes.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120323    D. Dirkx          Removed set functions; moved functionality to constructor.
 *
 *    References
 *
 *    Notes
 *
 */

#include <cmath>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <TudatCore/Mathematics/BasicMathematics/coordinateConversions.h>

#include "Tudat/Mathematics/GeometricShapes/sphereSegment.h"

namespace tudat
{
namespace geometric_shapes
{

// Using declarations.
using std::cerr;
using std::endl;
using std::cos;
using std::sin;
using tudat::unit_conversions::convertRadiansToDegrees;

SphereSegment::SphereSegment( const double radius,
                              const double minimumAzimuthAngle,
                              const double maximumAzimuthAngle,
                              const double minimumZenithAngle,
                              const double maximumZenithAngle )
{
    radius_ = radius;
    setMinimumIndependentVariable( 1, minimumAzimuthAngle );
    setMaximumIndependentVariable( 1, maximumAzimuthAngle );
    setMinimumIndependentVariable( 2, minimumZenithAngle );
    setMaximumIndependentVariable( 2, maximumZenithAngle );
}

//! Get surface point on sphere segment.
Eigen::VectorXd SphereSegment::getSurfacePoint( const double azimuthAngle,
                                                const double zenithAngle )
{
    Eigen::Vector3d sphericalPositionVector = Eigen::Vector3d( radius_, zenithAngle, azimuthAngle );

    // Gets surface point on sphere, unrotated and centered at origin.
    cartesianPositionVector_
            = basic_mathematics::coordinate_conversions::convertSphericalToCartesian(
                sphericalPositionVector );

    // Translate point.
    transformPoint( cartesianPositionVector_ );

    // Return Cartesian position vector.
    return cartesianPositionVector_;
}

//! Get surface derivative on sphere segment.
Eigen::VectorXd SphereSegment::getSurfaceDerivative( const double azimuthAngle,
                                                     const double zenithAngle,
                                                     const int powerOfAzimuthAngleDerivative,
                                                     const int powerOfZenithAngleDerivative )
{
    // Declare and set size of derivative vector.
    Eigen::VectorXd derivative_ = Eigen::VectorXd( 3 );

    // Go through the different possibilities for the values of the power
    // of the derivative.
    if ( powerOfAzimuthAngleDerivative < 0 || powerOfZenithAngleDerivative < 0 )
    {
        derivative_( 0 ) = 0.0;
        derivative_( 1 ) = 0.0;
        derivative_( 2 ) = 0.0;

        cerr << "No negative power of derivatives allowed, returning 0,0,0" << endl;
    }

    // When requesting the zeroth derivative with respect to the two
    // independent variables, the surface point is returned. Note that this
    // does include the offset.
    else if ( powerOfAzimuthAngleDerivative == 0 && powerOfZenithAngleDerivative == 0 )
    {
        derivative_ = getSurfacePoint( azimuthAngle, zenithAngle );
    }

    // The contributions of the two parameterizing variables are determined
    // and then multiplied component-wise.
    else
    {
        // Declare and set sizes contribution vectors.
        Eigen::VectorXd derivative1Contribution_ = Eigen::VectorXd( 3 );
        Eigen::VectorXd derivative2Contribution_ = Eigen::VectorXd( 3 );

        // Since this derivative is "cyclical", as it is
        // only dependent on sines and cosines, only the "modulo 4"th
        // derivatives need to be determined. Derivatives are determined
        // from the form of the spherical coordinates, see
        // basic_mathematics::coordinateConversions::convertSphericalToCartesian from
        // Tudat Core.
        switch( powerOfAzimuthAngleDerivative % 4 )
        {
        case( 0 ):

            derivative1Contribution_( 0 ) = cos( azimuthAngle );
            derivative1Contribution_( 1 ) = sin( azimuthAngle );
            derivative1Contribution_( 2 ) = 0.0;
            break;

        case( 1 ):

            derivative1Contribution_( 0 ) = -sin( azimuthAngle );
            derivative1Contribution_( 1 ) = cos( azimuthAngle );
            derivative1Contribution_( 2 ) = 0.0;
            break;

        case( 2 ):

            derivative1Contribution_( 0 ) = -cos( azimuthAngle );
            derivative1Contribution_( 1 ) = -sin( azimuthAngle );
            derivative1Contribution_( 2 ) = 0.0;
            break;

        case( 3 ):

            derivative1Contribution_( 0 ) = sin( azimuthAngle );
            derivative1Contribution_( 1 ) = -cos( azimuthAngle );
            derivative1Contribution_( 2 ) = 0.0;
            break;

        default:

            cerr << " Bad value for powerOfAzimuthAngleDerivative"
                 << " ( mod 4 ) of value is not 0, 1, 2 or 3 " << endl;
        }

        // This derivative is "cyclical" in the same manner as the derivative
        // with respect to the 1st independent variable.
        switch( powerOfZenithAngleDerivative %4 )
        {
        case( 0 ):

            derivative2Contribution_( 0 ) = sin( zenithAngle );
            derivative2Contribution_( 1 ) = sin( zenithAngle );
            derivative2Contribution_( 2 ) = cos( zenithAngle );
            break;

        case( 1 ):

            derivative2Contribution_( 0 ) = cos( zenithAngle );
            derivative2Contribution_( 1 ) = cos( zenithAngle );
            derivative2Contribution_( 2 ) = -sin( zenithAngle );
            break;

        case( 2 ):

            derivative2Contribution_( 0 ) = -sin( zenithAngle );
            derivative2Contribution_( 1 ) = -sin( zenithAngle );
            derivative2Contribution_( 2 ) = -cos( zenithAngle );
            break;

        case( 3 ):

            derivative2Contribution_( 0 ) = -cos( zenithAngle );
            derivative2Contribution_( 1 ) = -cos( zenithAngle );
            derivative2Contribution_( 2 ) = sin( zenithAngle );
            break;

        default:

            cerr << " Bad value for powerOfZenithAngleDerivative"
                 << " ( mod 4 ) of value is not 0, 1, 2 or 3 " << endl;
        }

        // Construct the full derivative.
        derivative_( 0 ) = derivative1Contribution_( 0 ) * derivative2Contribution_( 0 );
        derivative_( 1 ) = derivative1Contribution_( 1 ) * derivative2Contribution_( 1 );
        derivative_( 2 ) = derivative1Contribution_( 2 ) * derivative2Contribution_( 2 );

        // Scale vector by radius.
        derivative_ = derivative_ * radius_;

        // Rotate derivatives to take into account rotation of sphere.
        derivative_ = rotationMatrix_ * scalingMatrix_ * derivative_;
    }

    return derivative_;
}

//! Get parameter of sphere segment.
double SphereSegment::getParameter( const int index )
{
    // Check if parameter is radius.
    if ( index == 0 )
    {
        parameter_ = radius_;
    }

    // Else return cerr statement.
    else
    {
        parameter_ = -0.0;

        // Cerr statement.
        cerr << "Parameter does not exist" << endl;
    }

    // Return parameter.
    return parameter_;
}

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream, SphereSegment& sphereSegment )
{
    stream << "This is a sphere segment geometry." << endl;
    stream << "The range of the independent variables are: " << endl;
    stream << "Azimuth angle: "
           << convertRadiansToDegrees( sphereSegment.getMinimumIndependentVariable( 1 ) )
           << " degrees to "
           << convertRadiansToDegrees( sphereSegment.getMaximumIndependentVariable( 1 ) )
           << " degrees" << endl;
    stream << "Zenith angle: "
           << convertRadiansToDegrees( sphereSegment.getMinimumIndependentVariable( 2 ) )
           << " degrees to "
           << convertRadiansToDegrees( sphereSegment.getMaximumIndependentVariable( 2 ) )
           << " degrees" << endl;
    stream << "The radius is: " << sphereSegment.getRadius( ) << " meter." << endl;

    // Return stream.
    return stream;
}

} // namespace geometric_shapes
} // namespace tudat

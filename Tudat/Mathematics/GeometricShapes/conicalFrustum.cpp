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
 *      102511    D. Dirkx          First version of file.
 *      110120    D. Dirkx          Finalized for code check.
 *      110208    K. Kumar          Updated file header; corrected Doxygen comments; minor changes.
 *      110209    D. Dirkx          Minor changes.
 *      110209    K. Kumar          Minor changes.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120118    D. Gondelach      Implemented new convertCylindricalToCartesianCoordinates.
 *      120323    D. Dirkx          Removed set functions; moved functionality to constructor
 *
 *    References
 *
 *    Notes
 *
 */

#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Mathematics/GeometricShapes/conicalFrustum.h"

namespace tudat
{
namespace geometric_shapes
{

using tudat::basic_mathematics::mathematical_constants::PI;
using std::cerr;
using std::endl;
using std::sin;
using std::cos;

//! Conical frustum consructor, sets all shape parameters.

ConicalFrustum::ConicalFrustum( const double coneHalfAngle,
                                const double startRadius,
                                const double length,
                                const double minimumAzimuthAngle,
                                const double maximumAzimuthAngle )
{
    coneHalfAngle_ = coneHalfAngle;
    startRadius_= startRadius;
    length_ = length;
    setMinimumIndependentVariable( 1, minimumAzimuthAngle );
    setMaximumIndependentVariable( 1, maximumAzimuthAngle );
    setMinimumIndependentVariable( 2, 0 );
    setMaximumIndependentVariable( 2, 1 );
}

//! Get surface point on conical frustum.
Eigen::VectorXd ConicalFrustum::getSurfacePoint( double azimuthAngle, double lengthFraction )
{
    // Determines the radius of the cone at the given length fraction.
    double localRadius_ = startRadius_ + length_ * lengthFraction * std::tan( coneHalfAngle_ );


    // Set x and y coordinate of untransformed cone.
    cartesianPositionVector_ = basic_mathematics::coordinate_conversions::
            convertCylindricalToCartesian( Eigen::Vector3d( localRadius_, azimuthAngle, 0.0 ) );

    // Set z coordinate of untransformed cone.
    cartesianPositionVector_( 2 ) = -length_ * lengthFraction;

    // Transform conical frustum to desired position and orientation.
    transformPoint( cartesianPositionVector_ );

    // Return Cartesian position vector.
    return cartesianPositionVector_;
}

//! Get surface derivative on conical frustum.
Eigen::VectorXd ConicalFrustum::getSurfaceDerivative(
        const double lengthFraction, const double azimuthAngle,
        const int powerOfLengthFractionDerivative, const int powerOfAzimuthAngleDerivative )
{
    // Declare and set size of derivative vector.
    Eigen::VectorXd derivative_ = Eigen::VectorXd( 3 );

    // No negative derivatives may be retrieved, a zero vector is returned in
    // this case.
    if ( powerOfLengthFractionDerivative < 0 || powerOfAzimuthAngleDerivative < 0 )
    {
        derivative_( 0 ) = 0.0;
        derivative_( 1 ) = 0.0;
        derivative_( 2 ) = 0.0;

        cerr << " No negative derivatives allowed when retrieving cone "
             << "derivative, returning 0,0,0" << endl;
    }

    // When requesting the zeroth derivative with respect to the two
    // independent variables, the surface point is returned.
    else if ( powerOfLengthFractionDerivative == 0 && powerOfAzimuthAngleDerivative == 0 )
    {
        derivative_ = getSurfacePoint( lengthFraction, azimuthAngle );
    }

    // The contributions of the two parametrizing variables are determined and
    // then multiplied component-wise.
    else
    {
        // Declare and set sizes of contribution vectors.
        Eigen::VectorXd derivative1Contribution_ = Eigen::VectorXd( 3 );
        Eigen::VectorXd derivative2Contribution_ = Eigen::VectorXd( 3 );

        switch( powerOfLengthFractionDerivative )
        {
        case( 0 ):

            derivative1Contribution_( 0 ) = 0.0;
            derivative1Contribution_( 1 ) = 0.0;
            derivative1Contribution_( 2 ) = -length_ * lengthFraction;
            break;

        case( 1 ):

            derivative1Contribution_( 0 ) = 0.0;
            derivative1Contribution_( 1 ) = 0.0;
            derivative1Contribution_( 2 ) = -length_;
            break;

        // For all higher derivatives, all components are zero.
        default:

            derivative1Contribution_( 0 ) = 0.0;
            derivative1Contribution_( 1 ) = 0.0;
            derivative1Contribution_( 2 ) = 0.0;
            break;
        }

        // Since this derivative is "cyclical", as it is only dependant on sines
        // and cosines, only the "modulo 4"th derivative need be determined.
        // Derivatives are determined from the form of the cylindrical coordinates,
        // see basic_mathematics::coordinateConversions::convertSphericalToCartesian from
        // Tudat Core.
        switch( powerOfAzimuthAngleDerivative % 4 )
        {
        case( 0 ):

            derivative2Contribution_( 0 ) = -sin( azimuthAngle );
            derivative2Contribution_( 1 ) = cos( azimuthAngle );
            derivative2Contribution_( 2 ) = 0.0;
            break;

        case( 1 ):

            derivative2Contribution_( 0 ) = -cos( azimuthAngle );
            derivative2Contribution_( 1 ) = -sin( azimuthAngle );
            derivative2Contribution_( 2 ) = 0.0;
            break;

        case( 2 ):

            derivative2Contribution_( 0 ) = sin( azimuthAngle );
            derivative2Contribution_( 1 ) = -cos( azimuthAngle );
            derivative2Contribution_( 2 ) = 0.0;
            break;

        case( 3 ):

            derivative2Contribution_( 0 ) = cos( azimuthAngle );
            derivative2Contribution_( 1 ) = sin( azimuthAngle );
            derivative2Contribution_( 2 ) = 0.0;
            break;

        default:

            cerr << " Bad value for powerOfAzimuthAngleDerivative "
                 << " ( mod 4 ) of value is not 0, 1, 2 or 3 " << endl;
        }

        // Combine contributions to derivative.
        derivative_( 0 ) = derivative1Contribution_( 0 ) * derivative2Contribution_( 0 );
        derivative_( 1 ) = derivative1Contribution_( 1 ) * derivative2Contribution_( 1 );
        derivative_( 2 ) = derivative1Contribution_( 2 ) * derivative2Contribution_( 2 );

        // Rotate derivatives to take into account rotation of cone.
        derivative_ = scalingMatrix_ * rotationMatrix_ * derivative_;
    }

    // Return derivative vector.
    return derivative_;
}

//! Get parameter of conical frustum.
double ConicalFrustum::getParameter( const int index )
{
   // Check which parameter is selected.
   switch( index )
   {
   case( 0 ):

       parameter_ = startRadius_;
       break;

   case( 1 ):

       parameter_ = length_;
       break;

   case( 2 ):

       parameter_ = coneHalfAngle_;
       break;

   default:

       cerr << "Parameter " << index << " does not exist in ConicalFrustum, returning 0" << endl;
       parameter_ = 0;
       break;
   }

   // Return parameter.
   return parameter_;
}

//! Overload ostream to print class information.
std::ostream &operator<<( std::ostream &stream, ConicalFrustum& conicalFrustum )
{
    stream << "This is a conical frustum geometry." << endl;
    stream << "The circumferential angle runs from: "
           << conicalFrustum.getMinimumAzimuthAngle( ) * 180.0 / PI << " degrees to "
           << conicalFrustum.getMaximumAzimuthAngle( ) * 180.0 / PI << " degrees" << endl;
    stream << "The start radius is: " << conicalFrustum.getStartRadius( ) << endl;
    stream << "The length is: " << conicalFrustum.getLength( ) << endl;
    stream << "The cone half angle is: "  << conicalFrustum.getConeHalfAngle( ) * 180.0 / PI
           << " degrees" << endl;

    // Return stream.
    return stream;
}

} // namespace geometric_shapes
} // namespace tudat

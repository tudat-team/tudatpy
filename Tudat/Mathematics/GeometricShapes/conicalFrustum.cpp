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
 *      102511    D. Dirkx          First version of file.
 *      110120    D. Dirkx          Finalized for code check.
 *      110208    K. Kumar          Updated file header; corrected Doxygen comments; minor changes.
 *      110209    D. Dirkx          Minor changes.
 *      110209    K. Kumar          Minor changes.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120118    D. Gondelach      Implemented new convertCylindricalToCartesianCoordinates.
 *
 *    References
 *
 */

// Temporary notes (move to class/function doxygen):
// The getSurfacePoint currently uses a VectorXd as a return type,
// this could be changed to a CartesianPositionElements type in the
// future for consistency with the rest of the code.
// 

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Mathematics/GeometricShapes/conicalFrustum.h"

namespace tudat
{

// Using declarations.
using tudat::mathematics::PI;
using std::cerr;
using std::endl;
using std::sin;
using std::cos;

//! Get surface point on conical frustum.
Eigen::VectorXd ConicalFrustum::getSurfacePoint( double azimuthAngle, double lengthFraction )
{
    // Determines the radius of the cone at the given length fraction.
    double localRadius = startRadius_ + length_ * lengthFraction * std::tan( coneHalfAngle_ );

    // Compute x, y and z coordinates of untransformed cone.
    Eigen::Vector3d cartesianPosition = mathematics::coordinate_conversions::
            convertCylindricalToCartesian(
                localRadius, azimuthAngle, -length_ * lengthFraction );

    // Set Cartesian coordinates of untransformed cone.
    cartesianPositionVector_.head( 3 ) = cartesianPosition;

    // Transform conical frustum to desired position and orientation.
    transformPoint( cartesianPositionVector_ );

    // Return Cartesian position vector.
    return cartesianPositionVector_;
}

//! Get surface derivative on conical frustum.
Eigen::VectorXd ConicalFrustum::getSurfaceDerivative(
        double lengthFraction, double azimuthAngle,
        int powerOfLengthFractionDerivative, int powerOfAzimuthAngleDerivative )
{
    // Declare and set size of derivative vector.
    Eigen::VectorXd derivative_ = Eigen::VectorXd( 3 );

    // No negative derivatives may be retrieved, a zero vector is returned in
    // this case.
    if ( powerOfLengthFractionDerivative < 0 || powerOfAzimuthAngleDerivative < 0 )
    {
        derivative_( 0 ) = 0;
        derivative_( 1 ) = 0;
        derivative_( 2 ) = 0;

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

            derivative1Contribution_( 0 ) = 0;
            derivative1Contribution_( 1 ) = 0;
            derivative1Contribution_( 2 ) = -length_ * lengthFraction;
            break;

        case( 1 ):

            derivative1Contribution_( 0 ) = 0;
            derivative1Contribution_( 1 ) = 0;
            derivative1Contribution_( 2 ) = -length_;
            break;

        // For all higher derivatives, all components are zero.
        default:

            derivative1Contribution_( 0 ) = 0;
            derivative1Contribution_( 1 ) = 0;
            derivative1Contribution_( 2 ) = 0;
            break;
        }

        // Since this derivative is "cyclical", as it is only dependant on sines
        // and cosines, only the "modulo 4"th derivative need be determined.
        // Derivatives are determined from the form of the cylindrical
        // coordinates, see mathematics::convertSphericalToCartesian.
        switch( powerOfAzimuthAngleDerivative % 4 )
        {
        case( 0 ):

            derivative2Contribution_( 0 ) = -sin( azimuthAngle );
            derivative2Contribution_( 1 ) = cos( azimuthAngle );
            derivative2Contribution_( 2 ) = 0 ;
            break;

        case( 1 ):

            derivative2Contribution_( 0 ) = -cos( azimuthAngle );
            derivative2Contribution_( 1 ) = -sin( azimuthAngle );
            derivative2Contribution_( 2 ) = 0 ;
            break;

        case( 2 ):

            derivative2Contribution_( 0 ) = sin( azimuthAngle );
            derivative2Contribution_( 1 ) = -cos( azimuthAngle );
            derivative2Contribution_( 2 ) = 0;
            break;

        case( 3 ):

            derivative2Contribution_( 0 ) = cos( azimuthAngle );
            derivative2Contribution_( 1 ) = sin( azimuthAngle );
            derivative2Contribution_( 2 ) = 0;
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
double ConicalFrustum::getParameter( int index )
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

//! Set a parameter of conical frustum.
void ConicalFrustum::setParameter( int index, double parameter )
{
    // Check which parameter is to be set and set value.
    switch( index )
    {
    case( 0 ):

        startRadius_ = parameter;
        break;

    case( 1 ):

        length_ = parameter;
        break;

    case( 2 ):

        coneHalfAngle_ = parameter;
        break;

    default:

        cerr << "Parameter " << index << " does not exist in "
             << "ConicalFrustum, returning 0";
        break;
    }
}

//! Overload ostream to print class information.
std::ostream &operator<<( std::ostream &stream, ConicalFrustum& conicalFrustum )
{
    stream << "This is a conical frustum geometry." << endl;
    stream << "The circumferential angle runs from: "
           << conicalFrustum.getMinimumAzimuthAngle( ) * 180/PI << " degrees to "
           << conicalFrustum.getMaximumAzimuthAngle( ) * 180/PI << " degrees" << endl;
    stream << "The start radius is: " << conicalFrustum.getStartRadius( ) << endl;
    stream << "The length is: " << conicalFrustum.getLength( ) << endl;
    stream << "The cone half angle is: "  << conicalFrustum.getConeHalfAngle( ) * 180 / PI
           <<  " degrees" << endl;

    // Return stream.
    return stream;
}

} // namespace tudat

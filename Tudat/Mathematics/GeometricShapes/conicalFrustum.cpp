/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Mathematics/GeometricShapes/conicalFrustum.h"


namespace tudat
{
namespace geometric_shapes
{

using mathematical_constants::PI;
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
    cartesianPositionVector_ = coordinate_conversions::
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

        throw std::runtime_error( " No negative derivatives allowed when retrieving cone derivative " );
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
        // see basic_mathematics::coordinateConversions::convertSphericalToCartesian
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

            throw std::runtime_error( " Bad value for powerOfAzimuthAngleDerivative ( mod 4 ) of value is not 0, 1, 2 or 3 ");
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
       std::string errorMessage = "Parameter " + std::to_string( index ) + "does not exist in ConicalFrustum.";
       throw std::runtime_error( errorMessage );
   }

   // Return parameter.
   return parameter_;
}

//! Overload ostream to print class information.
std::ostream &operator << ( std::ostream &stream, ConicalFrustum& conicalFrustum )
{
    stream << "This is a conical frustum geometry." << std::endl;
    stream << "The circumferential angle runs from: "
           << conicalFrustum.getMinimumAzimuthAngle( ) * 180.0 / PI << " degrees to "
           << conicalFrustum.getMaximumAzimuthAngle( ) * 180.0 / PI << " degrees" << std::endl;
    stream << "The start radius is: " << conicalFrustum.getStartRadius( ) << std::endl;
    stream << "The length is: " << conicalFrustum.getLength( ) << std::endl;
    stream << "The cone half angle is: "  << conicalFrustum.getConeHalfAngle( ) * 180.0 / PI
           << " degrees" << std::endl;

    // Return stream.
    return stream;
}

} // namespace geometric_shapes
} // namespace tudat

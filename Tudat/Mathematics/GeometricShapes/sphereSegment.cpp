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

#include <cmath>

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"

#include "Tudat/Mathematics/GeometricShapes/sphereSegment.h"

namespace tudat
{
namespace geometric_shapes
{

// Using declarations.
using std::cos;
using std::sin;
using unit_conversions::convertRadiansToDegrees;

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
            = coordinate_conversions::convertSphericalToCartesian(
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

        throw std::runtime_error( "No negative power of derivatives allowed in sphere segment." );
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
        // basic_mathematics::coordinateConversions::convertSphericalToCartesian
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

            throw std::runtime_error( " Bad value for powerOfAzimuthAngleDerivative ( mod 4 ) of value is not 0, 1, 2 or 3 in sphere segment." );

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

            throw std::runtime_error( " Bad value for powerOfZenithAngleDerivative ( mod 4 ) of value is not 0, 1, 2 or 3 in sphere segment." );
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
        std::string errorMessage = "Parameter " + std::to_string( index ) + "does not exist in sphere segment.";
        throw std::runtime_error( errorMessage );
    }

    // Return parameter.
    return parameter_;
}

//! Overload ostream to print class information.
std::ostream& operator << ( std::ostream& stream, SphereSegment& sphereSegment )
{
    stream << "This is a sphere segment geometry." << std::endl;
    stream << "The range of the independent variables are: " << std::endl;
    stream << "Azimuth angle: "
           << convertRadiansToDegrees( sphereSegment.getMinimumIndependentVariable( 1 ) )
           << " degrees to "
           << convertRadiansToDegrees( sphereSegment.getMaximumIndependentVariable( 1 ) )
           << " degrees" << std::endl;
    stream << "Zenith angle: "
           << convertRadiansToDegrees( sphereSegment.getMinimumIndependentVariable( 2 ) )
           << " degrees to "
           << convertRadiansToDegrees( sphereSegment.getMaximumIndependentVariable( 2 ) )
           << " degrees" << std::endl;
    stream << "The radius is: " << sphereSegment.getRadius( ) << " meter." << std::endl;

    // Return stream.
    return stream;
}

} // namespace geometric_shapes
} // namespace tudat

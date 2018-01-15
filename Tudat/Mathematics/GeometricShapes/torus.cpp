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

#include <iostream>

#include "Tudat/Mathematics/GeometricShapes/torus.h"

namespace tudat
{
namespace geometric_shapes
{

using std::sin;
using std::cos;
using mathematical_constants::PI;

Torus::Torus( const double majorRadius, const double minorRadius,
              const double minimumMajorCircumferentialAngle,
              const double maximumMajorCircumferentialAngle,
              const double minimumMinorCircumferentialAngle,
              const double maximumMinorCircumferentialAngle )
{
    minorRadius_ = minorRadius;
    majorRadius_ = majorRadius;
            setMinimumIndependentVariable( 1, minimumMajorCircumferentialAngle );
            setMaximumIndependentVariable( 1, maximumMajorCircumferentialAngle );
            setMinimumIndependentVariable( 2, minimumMinorCircumferentialAngle );
            setMaximumIndependentVariable( 2, maximumMinorCircumferentialAngle );
}

//! Get surface point on torus.
Eigen::VectorXd Torus::getSurfacePoint( const double majorCircumferentialAngle,
                                        const double minorCircumferentialAngle )
{
    // Form cartesian position vector from paranmetrization.
    cartesianPositionVector_( 0 )
            = ( majorRadius_ + ( minorRadius_ * std::cos( minorCircumferentialAngle ) ) )
            * std::cos( majorCircumferentialAngle );
    cartesianPositionVector_( 1 )
            = ( majorRadius_ + ( minorRadius_ * std::cos( minorCircumferentialAngle ) ) )
            * std::sin( majorCircumferentialAngle );
    cartesianPositionVector_( 2 ) = minorRadius_ *  std::sin( minorCircumferentialAngle );

    // Translate vector.
    transformPoint( cartesianPositionVector_ );

    // Return Cartesian position vector.
    return cartesianPositionVector_;
}

//! Get surface derivative on torus.
Eigen::VectorXd Torus::getSurfaceDerivative( const double majorCircumferentialAngle,
                                             const double minorCircumferentialAngle,
                                             const int powerOfMajorCircumferentialAngleDerivative,
                                             const int powerOfMinorCircumferentialAngleDerivative )
{
    // Declare and set size of derivative vector.
    Eigen::VectorXd derivative_ = Eigen::VectorXd( 3 );

    // Go through the different possibilities for the values of the power
    // of the derivative.
    if ( powerOfMajorCircumferentialAngleDerivative < 0
         || powerOfMinorCircumferentialAngleDerivative < 0 )
    {
        derivative_( 0 ) = 0.0;
        derivative_( 1 ) = 0.0;
        derivative_( 2 ) = 0.0;
    }

    // When requesting the zeroth derivative with respect to the two
    // independent variables, the surface point is returned. Note that this
    // does include the offset.
    else if ( powerOfMajorCircumferentialAngleDerivative == 0 &&
              powerOfMinorCircumferentialAngleDerivative == 0 )
    {
        derivative_ = getSurfacePoint( majorCircumferentialAngle,
                                       minorCircumferentialAngle );
    }

    // The contributions of the two parameterizing variables are determined
    // and then multiplied component-wise.
    else
    {
        // Declare and set sizes of contribution vectors.
        Eigen::VectorXd derivative1Contribution_ = Eigen::VectorXd( 3 );
        Eigen::VectorXd derivative2Contribution_ = Eigen::VectorXd( 3 );

        // Since this derivative is "cyclical", as it is
        // only dependent on sines and cosines, only the "modulo 4"th
        // derivatives need to be determined. Derivatives are determined
        // from the form of the spherical coordinates, see
        // basic_mathematics::coordinateConversions::convertSphericalToCartesian
        switch( powerOfMajorCircumferentialAngleDerivative % 4 )
        {
        case( 0 ):

            derivative1Contribution_( 0 ) = std::cos( majorCircumferentialAngle );
            derivative1Contribution_( 1 ) = std::sin( majorCircumferentialAngle );
            derivative1Contribution_( 2 ) = 1.0;
            break;

        case( 1 ):

            derivative1Contribution_( 0 ) = -std::sin( majorCircumferentialAngle );
            derivative1Contribution_( 1 ) = std::cos( majorCircumferentialAngle );
            derivative1Contribution_( 2 ) = 0.0;
            break;

        case( 2 ):

            derivative1Contribution_( 0 ) = -std::cos( majorCircumferentialAngle );
            derivative1Contribution_( 1 ) = -std::sin( majorCircumferentialAngle );
            derivative1Contribution_( 2 ) = 0.0;
            break;

        case( 3 ):

            derivative1Contribution_( 0 ) = std::sin( majorCircumferentialAngle );
            derivative1Contribution_( 1 ) = -std::cos( majorCircumferentialAngle );
            derivative1Contribution_( 2 ) = 0.0;
            break;

        default:

            throw std::runtime_error( " Bad value for powerOfMajorCircumferentialAngleDerivative ( mod 4 ) of value is not 0, 1, 2 or 3 in torus" );
        }

        // This derivative is "cyclical" in the same manner as the derivative
        // with respect to the 1st independent variable. One check needs to
        // made, as the zeroth derivative violates cyclicity
        switch( powerOfMinorCircumferentialAngleDerivative %4 )
        {
        case( 0 ):

            derivative2Contribution_( 0 ) = minorRadius_ * std::cos( minorCircumferentialAngle );
            derivative2Contribution_( 1 ) = minorRadius_ * std::cos( minorCircumferentialAngle );
            derivative2Contribution_( 2 ) = minorRadius_ * std::sin( minorCircumferentialAngle );

            // For the zeroth derivative, a term needs to be added.
            if ( powerOfMinorCircumferentialAngleDerivative == 0 )
            {
                derivative2Contribution_( 0 ) += majorRadius_;
                derivative2Contribution_( 1 ) += majorRadius_;
            }

            break;

        case( 1 ):

            derivative2Contribution_( 0 ) = -minorRadius_ * std::sin( minorCircumferentialAngle );
            derivative2Contribution_( 1 ) = -minorRadius_ * std::sin( minorCircumferentialAngle );
            derivative2Contribution_( 2 ) = minorRadius_ * std::cos( minorCircumferentialAngle );
            break;

        case( 2 ):

            derivative2Contribution_( 0 ) = -minorRadius_ * std::cos( minorCircumferentialAngle );
            derivative2Contribution_( 1 ) = -minorRadius_ * std::cos( minorCircumferentialAngle );
            derivative2Contribution_( 2 ) = -minorRadius_ * std::sin( minorCircumferentialAngle );
            break;

        case( 3 ):

            derivative2Contribution_( 0 ) =  minorRadius_ * std::sin( minorCircumferentialAngle );
            derivative2Contribution_( 1 ) =  minorRadius_ * std::sin( minorCircumferentialAngle );
            derivative2Contribution_( 2 ) = -minorRadius_ * std::cos( minorCircumferentialAngle );
            break;

        default:

            throw std::runtime_error( " Bad value for powerOfMinorCircumferentialAngleDerivative ( mod 4 ) of value is not 0, 1, 2 or 3 in torus" );
        }

        // Construct the full derivative.
        derivative_( 0 ) = derivative1Contribution_( 0 ) * derivative2Contribution_( 0 );
        derivative_( 1 ) = derivative1Contribution_( 1 ) * derivative2Contribution_( 1 );
        derivative_( 2 ) = derivative1Contribution_( 2 ) * derivative2Contribution_( 2 );

        // Rotate derivatives to take into account rotation of torus.
        derivative_ = rotationMatrix_ * scalingMatrix_ * derivative_;
    }

    // Return derivative vector.
    return derivative_;
}

//! Get parameter of torus.
double Torus::getParameter( const int index )
{
    // Set parameter based on index.
    switch( index )
    {
    case 0:

        parameter_ = majorRadius_;
        break;

    case 1:

        parameter_ = minorRadius_;
        break;

    default:
        std::string errorMessage = "Parameter " + std::to_string( index ) + "does not exist in torus.";
        throw std::runtime_error( errorMessage );
    }

    // Return parameter.
    return parameter_;
}


//! Overload ostream to print class information.
std::ostream &operator << ( std::ostream &stream, Torus& torus )
{
    stream << "This is a torus geometry." << std::endl;
    stream << "The minor angle runs from: "
           << torus.getMinimumMinorCircumferentialAngle( ) * 180.0 / PI << " degrees to "
           << torus.getMaximumMinorCircumferentialAngle( ) * 180.0 / PI << " degrees" << std::endl;
    stream << "The major angle runs from: "
           << torus.getMinimumMajorCircumferentialAngle( ) * 180.0 / PI << " degrees to "
           << torus.getMaximumMajorCircumferentialAngle( ) * 180.0 / PI << " degrees" << std::endl;
    stream << "The major radius is: " << torus.getMajorRadius( ) << std::endl;
    stream << "The minor ( tube ) radius is: " << torus.getMinorRadius( ) << std::endl;

    // Return stream.
    return stream;
}

} // namespace geometric_shapes
} // namespace tudat

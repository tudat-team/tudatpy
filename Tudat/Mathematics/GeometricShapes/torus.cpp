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
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *
 *    References
 *
 */

// Temporary notes (move to class/function doxygen):
// The getSurfacePoint currently uses a VectorXd as a return type,
// this could be changed to a CartesianPositionElements type in the
// future for consistency with the rest of the code.
// 

#include <iostream>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>
#include "Tudat/Mathematics/GeometricShapes/torus.h"

namespace tudat
{

// Using declarations.
using tudat::mathematics::PI;
using std::cerr;
using std::endl;
using std::sin;
using std::cos;

//! Get surface point on torus.
Eigen::VectorXd Torus::getSurfacePoint( double majorCircumferentialAngle,
                                        double minorCircumferentialAngle )
{
    // Form cartesian position vector from paranmetrization.
    cartesianPositionVector_( 0 )
            = ( majorRadius_ + ( minorRadius_ * cos( minorCircumferentialAngle ) ) )
            * cos( majorCircumferentialAngle );
    cartesianPositionVector_( 1 )
            = ( majorRadius_ + ( minorRadius_ * cos( minorCircumferentialAngle ) ) )
            * sin( majorCircumferentialAngle );
    cartesianPositionVector_( 2 ) = minorRadius_ *  sin( minorCircumferentialAngle );

    // Translate vector.
    transformPoint( cartesianPositionVector_ );

    // Return Cartesian position vector.
    return cartesianPositionVector_;
}

//! Get surface derivative on torus.
Eigen::VectorXd Torus::getSurfaceDerivative( double majorCircumferentialAngle,
                                             double minorCircumferentialAngle,
                                             int powerOfMajorCircumferentialAngleDerivative,
                                             int powerOfMinorCircumferentialAngleDerivative )
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
        // mathematics::convertSphericalToCartesian.
        switch( powerOfMajorCircumferentialAngleDerivative % 4 )
        {
        case( 0 ):

            derivative1Contribution_( 0 ) = cos( majorCircumferentialAngle );
            derivative1Contribution_( 1 ) = sin( majorCircumferentialAngle );
            derivative1Contribution_( 2 ) = 1.0;
            break;

        case( 1 ):

            derivative1Contribution_( 0 ) = -sin( majorCircumferentialAngle );
            derivative1Contribution_( 1 ) = cos( majorCircumferentialAngle );
            derivative1Contribution_( 2 ) = 0.0;
            break;

        case( 2 ):

            derivative1Contribution_( 0 ) = -cos( majorCircumferentialAngle );
            derivative1Contribution_( 1 ) = -sin( majorCircumferentialAngle );
            derivative1Contribution_( 2 ) = 0.0;
            break;

        case( 3 ):

            derivative1Contribution_( 0 ) = sin( majorCircumferentialAngle );
            derivative1Contribution_( 1 ) = -cos( majorCircumferentialAngle );
            derivative1Contribution_( 2 ) = 0.0;
            break;

        default:

            cerr << " Bad value for powerOfMajorCircumferentialAngleDerivative"
                 << " ( mod 4 ) of value is not 0, 1, 2 or 3 " << endl;
        }

        // This derivative is "cyclical" in the same manner as the derivative
        // with respect to the 1st independent variable. One check needs to
        // made, as the zeroth derivative violates cyclicity
        switch( powerOfMinorCircumferentialAngleDerivative %4 )
        {
        case( 0 ):

            derivative2Contribution_( 0 ) = minorRadius_ * cos( minorCircumferentialAngle );
            derivative2Contribution_( 1 ) = minorRadius_ * cos( minorCircumferentialAngle );
            derivative2Contribution_( 2 ) = minorRadius_ * sin( minorCircumferentialAngle );

            // For the zeroth derivative, a term needs to be added.
            if ( powerOfMinorCircumferentialAngleDerivative == 0 )
            {
                derivative2Contribution_( 0 ) += majorRadius_;
                derivative2Contribution_( 1 ) += majorRadius_;
            }

            break;

        case( 1 ):

            derivative2Contribution_( 0 ) = -minorRadius_ * sin( minorCircumferentialAngle );
            derivative2Contribution_( 1 ) = -minorRadius_ * sin( minorCircumferentialAngle );
            derivative2Contribution_( 2 ) = minorRadius_ * cos( minorCircumferentialAngle );
            break;

        case( 2 ):

            derivative2Contribution_( 0 ) = -minorRadius_ * cos( minorCircumferentialAngle );
            derivative2Contribution_( 1 ) = -minorRadius_ * cos( minorCircumferentialAngle );
            derivative2Contribution_( 2 ) = -minorRadius_ * sin( minorCircumferentialAngle );
            break;

        case( 3 ):

            derivative2Contribution_( 0 ) =  minorRadius_ * sin( minorCircumferentialAngle );
            derivative2Contribution_( 1 ) =  minorRadius_ * sin( minorCircumferentialAngle );
            derivative2Contribution_( 2 ) = -minorRadius_ * cos( minorCircumferentialAngle );
            break;

        default:

            cerr << " Bad value for powerOfMinorCircumferentialAngleDerivative"
                 << " ( mod 4 ) of value is not 0, 1, 2 or 3 " << endl;
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
double Torus::getParameter( int index )
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

        cerr << "Parameter " << index << " does not exist in Torus, returning 0" << endl;
        break;
    }

    // Return parameter.
    return parameter_;
}

//! Set parameter of torus.
void Torus::setParameter( int index, double parameter)
{
    // Check which parameter is to be set and set value.
    switch( index )
    {
    case 0:

        majorRadius_ = parameter;
        break;

    case 1:

        minorRadius_ = parameter;
        break;

    default:

        cerr << "Parameter " << index << " does not exist in Torus.";
        break;
    }
}

//! Overload ostream to print class information.
std::ostream &operator<<( std::ostream &stream, Torus& torus )
{
    stream << "This is a torus geometry." << endl;
    stream << "The minor angle runs from: "
           << torus.getMinimumMinorCircumferentialAngle( ) * 180.0 / PI << " degrees to "
           << torus.getMaximumMinorCircumferentialAngle( ) * 180.0 / PI << " degrees" << endl;
    stream << "The major angle runs from: "
           << torus.getMinimumMajorCircumferentialAngle( ) * 180.0 / PI << " degrees to "
           << torus.getMaximumMajorCircumferentialAngle( ) * 180/PI << " degrees" << endl;
    stream << "The major radius is: " << torus.getMajorRadius( ) << endl;
    stream << "The minor ( tube ) radius is: " << torus.getMinorRadius( ) << endl;

    // Return stream.
    return stream;
}

} // namespace tudat

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
 *      102511    D. Dirkx          First version of file.
 *      110120    D. Dirkx          Finalized for code check.
 *      110208    K. Kumar          Updated file header; correct Doxygen comments; minor changes
 *                                  to functions.
 *      110209    D. Dirkx          Minor changes.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120323    D. Dirkx          Removed set functions; moved functionality to constructor
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *
 *    References
 *
 */

#include <iostream>

#include "Tudat/Mathematics/GeometricShapes/torus.h"

namespace tudat
{
namespace geometric_shapes
{

using std::cerr;
using std::endl;
using std::sin;
using std::cos;
using basic_mathematics::mathematical_constants::PI;

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
        // from Tudat Core.
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

        cerr << "Parameter " << index << " does not exist in Torus, returning 0" << endl;
        break;
    }

    // Return parameter.
    return parameter_;
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
           << torus.getMaximumMajorCircumferentialAngle( ) * 180.0 / PI << " degrees" << endl;
    stream << "The major radius is: " << torus.getMajorRadius( ) << endl;
    stream << "The minor ( tube ) radius is: " << torus.getMinorRadius( ) << endl;

    // Return stream.
    return stream;
}

} // namespace geometric_shapes
} // namespace tudat

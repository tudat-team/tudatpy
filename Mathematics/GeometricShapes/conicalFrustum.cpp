/*! \file conicalFrustum.cpp
 *    This file contains the implementation of the ConicalFrustum class.
 *
 *    Path              : /Mathematics/GeometricShapes/
 *    Version           : 5
 *    Check status      : Checked
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 25 November, 2010
 *    Last modified     : 9 February, 2011
 *
 *    References
 *
 *    Notes
 *      The getSurfacePoint currently uses a VectorXd as a return type,
 *      this could be changed to a CartesianPositionElements type in the
 *      future for consistency with the rest of the code.
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
 *      102511    D. Dirkx          First version of file.
 *      110120    D. Dirkx          Finalized for code check.
 *      110208    K. Kumar          Updated file header; corrected Doxygen
 *                                  comments; minor changes.
 *      110209    D. Dirkx          Minor changes.
 *      110209    K. Kumar          Minor changes.
 */

// Include statements.
#include "conicalFrustum.h"

// Using declarations.
using std::cerr;
using std::endl;

//! Default constructor.
ConicalFrustum::ConicalFrustum( )
{
    minimumIndependentVariable2_ = 0;
    maximumIndependentVariable2_ = 1;
}

//! Default destructor.
ConicalFrustum::~ConicalFrustum( )
{

}

//! Get surface point on conical frustum.
VectorXd ConicalFrustum::getSurfacePoint( const double& azimuthAngle,
                                          const double& lengthFraction )
{
    // Determines the radius of the cone at the given length fraction.
    double localRadius_ = startRadius_ + length_ * lengthFraction *
                          tan( coneHalfAngle_ );

    // Set x and y coordinate of untransformed cone.
    mathematics::convertCylindricalToCartesian( localRadius_,
                                                azimuthAngle,
                                                cartesianPositionVector_ );

    // Set z coordinate of untransformed cone.
    cartesianPositionVector_( 2 ) = -length_ * lengthFraction;

    // Transform conical frustum to desired position and orientation.
    transformPoint( cartesianPositionVector_ );

    // Return Cartesian position vector.
    return cartesianPositionVector_;
}

//! Get surface derivative on conical frustum.
VectorXd ConicalFrustum::getSurfaceDerivative(
        const double& lengthFraction, const double& azimuthAngle,
        const int& powerOfLengthFractionDerivative,
        const int& powerOfAzimuthAngleDerivative )
{
    // Declare and set size of derivative vector.
    VectorXd derivative_ = VectorXd( 3 );

    // No negative derivatives may be retrieved, a zero vector is returned in
    // this case.
    if ( powerOfLengthFractionDerivative < 0 ||
         powerOfAzimuthAngleDerivative < 0 )
    {
        derivative_( 0 ) = 0;
        derivative_( 1 ) = 0;
        derivative_( 2 ) = 0;

        cerr << " No negative derivatives allowed when retrieving cone "
             << "derivative, returning 0,0,0" << endl;
    }

    // When requesting the zeroth derivative with respect to the two
    // independent variables, the surface point is returned.
    else if ( powerOfLengthFractionDerivative == 0 &&
              powerOfAzimuthAngleDerivative == 0 )
    {
        derivative_ = getSurfacePoint( lengthFraction, azimuthAngle );
    }

    // The contributions of the two parametrizing variables are determined and
    // then multiplied component-wise.
    else
    {
        // Declare and set sizes of contribution vectors.
        VectorXd derivative1Contribution_ = VectorXd( 3 );
        VectorXd derivative2Contribution_ = VectorXd( 3 );

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
        derivative_( 0 ) = derivative1Contribution_( 0 ) *
                           derivative2Contribution_( 0 );
        derivative_( 1 ) = derivative1Contribution_( 1 ) *
                           derivative2Contribution_( 1 );
        derivative_( 2 ) = derivative1Contribution_( 2 ) *
                           derivative2Contribution_( 2 );

        // Rotate derivatives to take into account rotation of cone.
        derivative_ = scalingMatrix_ * rotationMatrix_ * derivative_;
    }

    // Return derivative vector.
    return derivative_;
}

//! Get parameter of conical frustum.
double ConicalFrustum::getParameter( const int& index )
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

       cerr << "Parameter " << index << " does not exist in  "
            << "ConicalFrustum, returning 0" << endl;
       parameter_ = 0;
       break;
   }

   // Return parameter.
   return parameter_;
}

//! Set a parameter of conical frustum.
void ConicalFrustum::setParameter( const int& index, const double& parameter )
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

//! Set cone half angle.
void ConicalFrustum::setConeHalfAngle( const double& coneHalfAngle )
{
    coneHalfAngle_ = coneHalfAngle;
}
//! Set length.
//! Function to set the length.
void ConicalFrustum::setLength( const double& length )
{
    length_ = length;
}

//! Set start radius.
void ConicalFrustum::setStartRadius( const double& startRadius )
{
    startRadius_ = startRadius;
}

//! Set minimum azimuth angle.
void ConicalFrustum::setMinimumAzimuthAngle( const double&
                                             minimumAzimuthAngle )
{
    setMinimumIndependentVariable( 1, minimumAzimuthAngle );
}

//! Set maximum azimuth angle.
void ConicalFrustum::setMaximumAzimuthAngle( const double&
                                             maximumAzimuthAngle )
{
    setMaximumIndependentVariable( 1, maximumAzimuthAngle );
}

//! Get cone half angle.
double& ConicalFrustum::getConeHalfAngle( )
{
    return coneHalfAngle_;
}

double& ConicalFrustum::getLength( )
{
    return length_;
}

//! Get start radius.
double& ConicalFrustum::getStartRadius( )
{
    return startRadius_;
}

//! Get minimum azimuth angle.
double ConicalFrustum::getMinimumAzimuthAngle( )
{
    return minimumIndependentVariable1_;
}

//! Get maximum azimuth angle.
double ConicalFrustum::getMaximumAzimuthAngle( )
{
    return maximumIndependentVariable1_;
}

//! Overload ostream to print class information.
std::ostream &operator<<( std::ostream &stream,
                          ConicalFrustum& conicalFrustum )
{
    stream << "This is a conical frustum geometry." << endl;
    stream << "The circumferential angle runs from: "
           << conicalFrustum.getMinimumAzimuthAngle( ) * 180/M_PI
           << " degrees to "
           << conicalFrustum.getMaximumAzimuthAngle( ) * 180/M_PI
           << " degrees" << endl;
    stream << "The start radius is: " << conicalFrustum.getStartRadius( )
            << endl;
    stream << "The length is: " << conicalFrustum.getLength( ) << endl;
    stream << "The cone half angle is: "
            << conicalFrustum.getConeHalfAngle( ) * 180 / M_PI
            <<  " degrees" << endl;

    // Return stream.
    return stream;
}

// End of file.

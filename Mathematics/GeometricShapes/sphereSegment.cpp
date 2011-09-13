/*! \file sphereSegment.cpp
 *    This file contains the implementation of the Sphere Segment class.
 *
 *    Path              : /Mathematics/GeometricShapes/
 *    Version           : 8
 *    Check status      : Checked
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 29 September, 2010
 *    Last modified     : 9 February, 2011
 *
 *    References
 *
 *    Notes
 *      The getSurfacePoint currently uses a VectorXd as a return type,
 *      this could be changed to a CartesianPositionElements type in the
 *      future for consistency with the rest of the code.
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 */

// Include statements.
#include "sphereSegment.h"

// Using declarations.
using std::cerr;
using std::endl;

//! Default constructor.
SphereSegment::SphereSegment( ) : radius_( -0.0 )
{
}

//! Default destructor.
SphereSegment::~SphereSegment( )
{
}

//! Get surface point on sphere segment.
VectorXd SphereSegment::getSurfacePoint( const double& azimuthAngle,
                                         const double& zenithAngle )
{
    // Gets surface point on sphere, unrotated and centered at origin.
    mathematics::convertSphericalToCartesian( radius_, azimuthAngle,
                                              zenithAngle,
                                              cartesianPositionVector_ );

    // Translate point.
    transformPoint( cartesianPositionVector_ );

    // Return Cartesian position vector.
    return cartesianPositionVector_;
}

//! Get surface derivative on sphere segment.
VectorXd SphereSegment::getSurfaceDerivative(
        const double& azimuthAngle, const double& zenithAngle,
        const int& powerOfAzimuthAngleDerivative,
        const int& powerOfZenithAngleDerivative )
{
    // Declare and set size of derivative vector.
    VectorXd derivative_ = VectorXd( 3 );

    // Go through the different possibilities for the values of the power
    // of the derivative.
    if ( powerOfAzimuthAngleDerivative < 0
         || powerOfZenithAngleDerivative < 0 )
    {
        derivative_( 0 ) = 0.0;
        derivative_( 1 ) = 0.0;
        derivative_( 2 ) = 0.0;

        cerr << "No negative power of derivatives allowed, "
             << "returning 0,0,0" << endl;
    }

    // When requesting the zeroth derivative with respect to the two
    // independent variables, the surface point is returned. Note that this
    // does include the offset.
    else if ( powerOfAzimuthAngleDerivative == 0
              && powerOfZenithAngleDerivative == 0 )
    {
        derivative_ = getSurfacePoint( azimuthAngle, zenithAngle );
    }

    // The contributions of the two parameterizing variables are determined
    // and then multiplied component-wise.
    else
    {
        // Declare and set sizes contribution vectors.
        VectorXd derivative1Contribution_ = VectorXd( 3 );
        VectorXd derivative2Contribution_ = VectorXd( 3 );

        // Since this derivative is "cyclical", as it is
        // only dependent on sines and cosines, only the "modulo 4"th
        // derivatives need to be determined. Derivatives are determined
        // from the form of the spherical coordinates, see
        // mathematics::convertSphericalToCartesian.
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
        derivative_( 0 ) = derivative1Contribution_( 0 )
                           * derivative2Contribution_( 0 );
        derivative_( 1 ) = derivative1Contribution_( 1 )
                           * derivative2Contribution_( 1 );
        derivative_( 2 ) = derivative1Contribution_( 2 )
                           * derivative2Contribution_( 2 );

        // Scale vector by radius.
        derivative_ = derivative_ * radius_;

        // Rotate derivatives to take into account rotation of sphere.
        derivative_ = rotationMatrix_ * scalingMatrix_ * derivative_;
    }

    return derivative_;
}

//! Get parameter of sphere segment.
double SphereSegment::getParameter( const int& index )
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

//! Set parameter of sphere segment.
void SphereSegment::setParameter( const int& index, const double& parameter )
{
    if ( index == 0 )
    {
        radius_ = parameter;
    }
}

//! Get radius.
double& SphereSegment::getRadius( )
{
    return radius_;
}

//! Set radius.
void SphereSegment::setRadius( const double& radius )
{
    radius_ = radius;
}

//! Set maximum value of azimuth angle.
void SphereSegment::setMaximumAzimuthAngle( const double& maximumAzimuthAngle )
{
    setMaximumIndependentVariable( 1, maximumAzimuthAngle );
}

//! Set minimum value of azimuth angle.
void SphereSegment::setMinimumAzimuthAngle( const double& minimumAzimuthAngle )
{
    setMinimumIndependentVariable( 1, minimumAzimuthAngle );
}

//! Set maximum value of zenith angle.
void SphereSegment::setMaximumZenithAngle( const double& maximumZenithAngle )
{
    setMaximumIndependentVariable( 2, maximumZenithAngle );
}

//! Set minimum value of zenith angle.
void SphereSegment::setMinimumZenithAngle( const double& minimumZenithAngle )
{
    setMinimumIndependentVariable( 2, minimumZenithAngle );
}

//! Get maximum value of azimuth angle.
double SphereSegment::getMaximumAzimuthAngle( )
{
    return maximumIndependentVariable1_;
}

//! Get minimum value of azimuth angle.
double SphereSegment::getMinimumAzimuthAngle( )
{
    return minimumIndependentVariable1_;

}

//! Get maximum value of zenith angle.
double SphereSegment::getMaximumZenithAngle( )
{
    return maximumIndependentVariable2_;
}

//! Get minimum value of the zenith angle.
double SphereSegment::getMinimumZenithAngle( )
{
    return minimumIndependentVariable2_;
}

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream, SphereSegment& sphereSegment )
{
    stream << "This is a sphere segment geometry." << endl;
    stream << "The range of the independent variables are: " << endl;
    stream << "Azimuth angle: "
           << convertRadiansToDegrees( sphereSegment
                                       .getMinimumIndependentVariable( 1 ) )
           << " degrees to "
           << convertRadiansToDegrees( sphereSegment
                                       .getMaximumIndependentVariable( 1 ) )
           << " degrees" << endl;
    stream << "Zenith angle: "
           << convertRadiansToDegrees( sphereSegment
                                       .getMinimumIndependentVariable( 2 ) )
           << " degrees to "
           << convertRadiansToDegrees( sphereSegment
                                       .getMaximumIndependentVariable( 2 ) )
           << " degrees" << endl;
    stream << "The radius is: " << sphereSegment.getRadius( ) << " meter."
           << endl;

    // Return stream.
    return stream;
}

// End of file.

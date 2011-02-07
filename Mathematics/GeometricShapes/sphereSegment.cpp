/*! \file sphereSegment.cpp
 *  This file contains the definition of the Sphere Segment class
 *
 *  Path              : Mathematics/Geometry/
 *  Version           : 1
 *  Check status      : Checked
 *
 *  Author            : Dominic Dirkx
 *  Affiliation       : TU Delft
 *  E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *  Checker           : J. Melman
 *  Affiliation       : Delft University of Technology
 *  E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *  Date created      : 29 September, 2010
 *  Last modified     : 29 September, 2010
 *
 *  References
 *
 *  Notes
 *    The setParameter and getParameter functions could benefit from
 *    taking an enum as input instead of an int.
 *
 *    The original file was split into several parts, so that each class is
 *    defined separately in a file. Hence, the this file was created after the
 *    first entries in the changelog and the file version is lower than the
 *    number of entries in changelog.
 *
 *  Copyright (c) 2010 Delft University of Technology.
 *
 *  This software is protected by national and international copyright.
 *  Any unauthorized use, reproduction or modification is unlawful and
 *  will be prosecuted. Commercial and non-private application of the
 *  software in any form is strictly prohibited unless otherwise granted
 *  by the authors.
 *
 *  The code is provided without any warranty; without even the implied
 *  warranty of merchantibility or fitness for a particular purpose.
 *
 *  Changelog
 *    100910   D. Dirkx                    First version of file
 *    100915   D. Dirkx                    Modified to correct comments, 80
 *                                         lines rule, etc.
 *    100928   D. Dirkx                    Modifications following first
 *                                         checking iteration.
 *    100929   D. Dirkx                    Creation of separate file for class
 *
 */

#include "sphereSegment.h"

//! Default constructor.
SphereSegment::SphereSegment( )
{
}

//! Default destructor.
SphereSegment::~SphereSegment( )
{
}

//! Constructor to initialize parameters.
SphereSegment::SphereSegment( const double& minimumIndependentVariable1,
                              const double& maximumIndependentVariable1,
                              const double& minimumIndependentVariable2,
                              const double& maximumIndependentVariable2,
                              const double& radius)
{
    minimumIndependentVariable1_ = minimumIndependentVariable1;
    maximumIndependentVariable1_ = maximumIndependentVariable1;
    minimumIndependentVariable2_ = minimumIndependentVariable2;
    maximumIndependentVariable2_ = maximumIndependentVariable2;
    radius_ = radius;
}

 //! Function to retrieve a surface point on the sphere.
VectorXd SphereSegment::getSurfacePoint( const double& independentVariable1,
                                         const double& independentVariable2 )
{
    VectorXd cartesianPosition = VectorXd( 3 );

    // Gets surface point on sphere, unrotated and centered at origin.
    mathematics::convertSphericalToCartesian( radius_, independentVariable1,
                                              independentVariable2,
                                              cartesianPosition );
    VectorXd point = VectorXd( 3 );

    // Rotates point to correct orientation.
    point( 0 ) = rotationMatrix_( 0, 0 ) * cartesianPosition( 0 ) +
                 rotationMatrix_( 0, 1 ) * cartesianPosition( 1 ) +
                 rotationMatrix_( 0, 2 ) * cartesianPosition( 2 );
    point( 1 ) = rotationMatrix_( 1, 0 ) * cartesianPosition( 0 ) +
                 rotationMatrix_( 1, 1 ) * cartesianPosition( 1 ) +
                 rotationMatrix_( 1, 2 ) * cartesianPosition( 2 );
    point( 2 ) = rotationMatrix_( 2, 0 ) * cartesianPosition( 0 ) +
                 rotationMatrix_( 2, 1 ) * cartesianPosition( 1 ) +
                 rotationMatrix_( 2, 2 ) * cartesianPosition( 2 );

    // Translates point such that sphere origin is correctly located.
    point = point + offset_;

    return point;
}

//! Calculates the surface derivative on the sphere w.r.t. the independent
//! variables.
VectorXd SphereSegment::getSurfaceDerivative( const double& independentVariable1,
                                            const double& independentVariable2,
                                            const int& powerOfDerivative1,
                                            const int& powerOfDerivative2 )
{
    // Set size of derivative.
    VectorXd derivative = VectorXd( 3 );

    // Go through the different possibilities for the values of the power of the derivative.
    if ( powerOfDerivative1 < 0 || powerOfDerivative2 < 0 )
    {
        std::cout << "No negative power of derivatives allowed, returning 0,0,0" << std::endl;
        derivative( 0 ) = 0.0;
        derivative( 1 ) = 0.0;
        derivative( 2 ) = 0.0;
    }

    // When requesting the zeroth derivative w.r.t. the two independent variables,
    // the surface point is returned. Note that this does include the offset.
    else if ( powerOfDerivative1 == 0 && powerOfDerivative2 == 0 )
    {
        derivative = getSurfacePoint( independentVariable1, independentVariable2 );
    }

    // The contributions of the two parameterizing variables are determined
    // and then  component-wise multiplied.
    else
    {
        VectorXd derivative1Contribution = VectorXd( 3 );
        VectorXd derivative2Contribution = VectorXd( 3 );

        // Since this derivative is "cyclical", as it is
        // only dependent on sines and cosines, only the "modulo 4"th
        // derivatives need to be determined. Derivatives are determined
        // from the form of the spherical coordinates, see
        // mathematics::convertSphericalToCartesian.
        switch( powerOfDerivative1 % 4 )
        {
        case( 0 ):
            derivative1Contribution( 0 ) = cos( independentVariable1 );
            derivative1Contribution( 1 ) = sin( independentVariable1 );
            derivative1Contribution( 2 ) = 0.0;
            break;
        case( 1 ):
            derivative1Contribution( 0 ) = -sin( independentVariable1 );
            derivative1Contribution( 1 ) = cos( independentVariable1 );
            derivative1Contribution( 2 ) = 0.0;
            break;
        case( 2 ):
            derivative1Contribution( 0 ) = -cos( independentVariable1 );
            derivative1Contribution( 1 ) = -sin( independentVariable1 );
            derivative1Contribution( 2 ) = 0.0;
            break;
        case( 3 ):
            derivative1Contribution( 0 ) = sin( independentVariable1 );
            derivative1Contribution( 1 ) = -cos( independentVariable1 );
            derivative1Contribution( 2 ) = 0.0;
            break;
        }

        // This derivative is "cyclical" in the same manner as the derivative
        // w.r.t. the 1st independent variable.
        switch( powerOfDerivative2 %4)
        {
        case( 0 ):
            derivative2Contribution( 0 ) = sin( independentVariable2 );
            derivative2Contribution( 1 ) = sin( independentVariable2 );
            derivative2Contribution( 2 ) = cos( independentVariable2 );
            break;
        case( 1 ):
            derivative2Contribution( 0 ) = cos( independentVariable2 );
            derivative2Contribution( 1 ) = cos( independentVariable2 );
            derivative2Contribution( 2 ) = -sin( independentVariable2 );
            break;
        case( 2 ):
            derivative2Contribution( 0 ) = -sin( independentVariable2 );
            derivative2Contribution( 1 ) = -sin( independentVariable2 );
            derivative2Contribution( 2 ) = -cos( independentVariable2 );
            break;
        case( 3 ):
            derivative2Contribution( 0 ) = -cos( independentVariable2 );
            derivative2Contribution( 1 ) = -cos( independentVariable2 );
            derivative2Contribution( 2 ) = sin( independentVariable2 );
            break;
        }

        // Construct the full derivative.
        derivative( 0 ) = derivative1Contribution( 0 ) * derivative2Contribution( 0 );
        derivative( 1 ) = derivative1Contribution( 1 ) * derivative2Contribution( 1 );
        derivative( 2 ) = derivative1Contribution( 2 ) * derivative2Contribution( 2 );

        // Scale vector by radius
        derivative = derivative * radius_;

        // Rotate derivatives to take into account rotation of sphere
        VectorXd newDerivative = VectorXd( 3 );
        newDerivative( 0 ) = rotationMatrix_( 0, 0 ) * derivative( 0 ) +
                             rotationMatrix_( 0, 1 ) * derivative( 1 ) +
                             rotationMatrix_( 0, 2 ) * derivative( 2 );
        newDerivative( 1 ) = rotationMatrix_( 1, 0 ) * derivative( 0 ) +
                             rotationMatrix_( 1, 1 ) * derivative( 1 ) +
                             rotationMatrix_( 1, 2 ) * derivative( 2 );
        newDerivative( 2 ) = rotationMatrix_( 2, 0 ) * derivative( 0 ) +
                             rotationMatrix_( 2, 1 ) * derivative( 1 ) +
                             rotationMatrix_( 2, 2 ) * derivative( 2 );
        derivative = newDerivative;
    }

    return derivative;
}

//! Function to retrieve a parameter value (radius).
double SphereSegment::getParameter( const int& parameterIndex )
{
    double parameter;

    // For parameterIndex = 0 the radius is returned.
    if( parameterIndex == 0 )
    {
        parameter = radius_;
    }

    // For all other values of parameterIndex, no parameter exists, 0 is returned and a
    // warning is displayed.
    else
    {
        std::cout<< "Parameter with index " << parameterIndex << " does not exist in SphereSegment; returning 0.0.";
        parameter = 0.0;
    }

    return parameter;
}

//! Function to set a parameter's value (radius).
void SphereSegment::setParameter( const int & parameterIndex, const double & parameter)
{
    // For parameterIndex = 0 the radius is set.
    if( parameterIndex == 0 )
    {
        radius_ = parameter;
    }

    // For all other values of parameterIndex, no parameter exists and a warning is displayed.
    else
    {
        std::cout << "Parameter " << parameterIndex << " does not exist in SphereSegment.";
    }
}

//! Overloaded ostream to print class information.
std::ostream& operator<<( std::ostream& stream, SphereSegment& geometry )
{
    stream << "This is a sphere segment geometry." << std::endl;
    stream << "The range of the independent variables are: " << std::endl;
    stream << "Azimuth angle: " << convertRadiansToDegrees( geometry.getMinimumIndependentVariable( 1 ) )
           << " degrees to " << convertRadiansToDegrees( geometry.getMaximumIndependentVariable( 1 ) ) <<
              " degrees" << std::endl;
    stream << "Zenith angle: " << convertRadiansToDegrees( geometry.getMinimumIndependentVariable( 2 ) ) <<
              " degrees to " << convertRadiansToDegrees( geometry.getMaximumIndependentVariable( 2 ) ) <<
              " degrees" << std::endl;
    stream << "The radius is: " << geometry.getParameter( 0 ) << " meter." << std::endl;

    return stream;
}

// End of file

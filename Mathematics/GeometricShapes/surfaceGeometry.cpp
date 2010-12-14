/*! \file surfaceGeometry.h
 *  This file contains the definition of the Surface Geometry base class
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

#include "surfaceGeometry.h"

//! Default constructor.
SurfaceGeometry::SurfaceGeometry( )
{
    // Declare the size of the offset vector.
    offset_.setZero( 3 );

    // Declare the rotation matrix.
    rotationMatrix_.setIdentity( 3, 3 );
}

//! Default destructor.
SurfaceGeometry::~SurfaceGeometry( )
{
}

//! Function to set the rotation and translation values.
void SurfaceGeometry::setTransFormValues( const VectorXd& offset,
                                          const MatrixXd& rotationMatrix )
{
    offset_ = offset;
    rotationMatrix_ = rotationMatrix;
}

//! Function to retrieve rotationMatrix_ from the shape.
MatrixXd SurfaceGeometry::getRotationMatrix( )
{
    return rotationMatrix_;
}

//! Function to retrieve offset from the shape.
VectorXd SurfaceGeometry::getOffset( )
{
    return offset_;
}

//! Function to set maximum values of independent variables.
void SurfaceGeometry::setMaximumIndependentVariable( const int& parameterIndex,
                                                     const double& value )
{
    if( parameterIndex == 1 )
    {
        maximumIndependentVariable1_ = value;
    }
    else if( parameterIndex == 2 )
    {
        maximumIndependentVariable2_ = value;
    }
    else
    {
        std::cout << "Parameter index " << parameterIndex << " not available in"
                  << "surface geometry; no value has been set.";
    }
}

//! Function to set minimum values of independent variables.
void SurfaceGeometry::setMinimumIndependentVariable( const int& parameterIndex,
                                                     const double& value )
{
    if( parameterIndex == 1 )
    {
        minimumIndependentVariable1_ = value;
    }
    else if( parameterIndex == 2 )
    {
        minimumIndependentVariable2_ = value;
    }
    else
    {
        std::cout << "Parameter index " << parameterIndex << " not available in"
                  << "surface geometry; no value has been set.";
    }
}

//! Function to get minimum values of independent variables.
double SurfaceGeometry::getMinimumIndependentVariable( const int& parameterIndex )
{
    double minimumValue;

    if( parameterIndex == 1 )
    {
        minimumValue = minimumIndependentVariable1_;
    }
    else if( parameterIndex == 2 )
    {
        minimumValue = minimumIndependentVariable2_;
    }
    else
    {
        std::cout << "Parameter index " << parameterIndex << " not available in"
                  << "surface geometry; returning 0.";
        minimumValue = 0.0;
    }
    return minimumValue;
}


//! Function to get maximum values of independent variables.
double SurfaceGeometry::getMaximumIndependentVariable( const int& parameterIndex )
{
    double maximumValue;

    if( parameterIndex == 1 )
    {
        maximumValue = maximumIndependentVariable1_;
    }
    else if( parameterIndex == 2 )
    {
        maximumValue = maximumIndependentVariable2_;
    }
    else
    {
        std::cout << "Parameter index " << parameterIndex << " not available in"
                  << "surface geometry; returning 0.";
        maximumValue = 0.0;
    }
    return maximumValue;
}

// End of file

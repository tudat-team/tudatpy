/*! \file singleSurfaceGeometry.cpp
 *    This file contains the definition of the SingleSurfaceGeometry base class.
 *
 *    Path              : /Mathematics/GeometricShapes/
 *    Version           : 10
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
 *    Last modified     : 5 September, 2011
 *
 *    References
 *
 *    Notes
 *      Contents of this file used to be in singleGeometry.cpp, but as this class
 *      has been split into single and composite surface geometry, the contents
 *      have been moved, with most of the SurfaceGeometry class now belonging to
 *      the SingleSurfaceGeometry class.
 *
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
 *      101125    D. Dirkx          Migration of contents to this file.
 *      110105    D. Dirkx          Final modifications for code check.
 *      110124    K. Kumar          Minor comment and layout changes.
 *      110204    K. Kumar          Minor comment and layout modifications;
 *                                  corrected Doxygen comments.
 *      110209    D. Dirkx          Minor changes.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

// Include statements.
#include "Mathematics/GeometricShapes/singleSurfaceGeometry.h"

//! Tudat library namespace.
namespace tudat
{

// Using declarations.
using std::cerr;
using std::endl;

//! Set minimum value of independent variable.
void SingleSurfaceGeometry::setMinimumIndependentVariable( const int& parameterIndex,
                                                           const double& minimumValue )
{
    independentVariable_ = IndependentVariables( parameterIndex );

    // Check which variable is to be set.
    switch ( independentVariable_ )
    {
    case firstIndependentVariable:

        // Set value.
        minimumIndependentVariable1_ = minimumValue;

        break;

    case secondIndependentVariable:

        // Set value.
        minimumIndependentVariable2_ = minimumValue;

        break;

    default:

        cerr << " Only 2 independent variables, variable "
             << parameterIndex << " does not exist when"
             << "setting minimum value"<< endl;
    }
}

//! Set maximum value of independent variable.
void SingleSurfaceGeometry::setMaximumIndependentVariable(
    const int& parameterIndex, const double& maximumValue )
{

    independentVariable_ = IndependentVariables( parameterIndex );

    // Check which variable is to be set.
    switch( independentVariable_ )
    {
    case firstIndependentVariable:

        // Set value.
        maximumIndependentVariable1_ = maximumValue;

        break;

    case secondIndependentVariable:

        // Set value.
        maximumIndependentVariable2_ = maximumValue;

        break;

    default:

        cerr << " Only 2 independent variables, variable "
             << parameterIndex << " does not exist when"
             << "setting maximum value"<< endl;
    }
}

//! Get minimum value of independent variable.
double SingleSurfaceGeometry::getMinimumIndependentVariable(
    const int& parameterIndex )
{
    // Declare local variables.
    double minimumValue_;

    // Check which variable is to be returned.
    switch( parameterIndex )
    {
    case 1:

        minimumValue_ = minimumIndependentVariable1_;

        break;

    case 2:

        minimumValue_ = minimumIndependentVariable2_;

        break;

    default:

        minimumValue_ = -0.0;

        cerr << " Only 2 independent variables, variable "
             <<  parameterIndex << " does not exist when "
             << "getting minimum value, returning -0.0" << endl;

    }

    // Return minimum value.
    return minimumValue_;
}

//! Get maximum value of independent variable.
double SingleSurfaceGeometry::getMaximumIndependentVariable(
    const int& parameterIndex )
{
    // Declare local variables.
    double maximumValue_;

    // Check which variable is to be returned.
    switch( parameterIndex )
    {
    case 1:

        maximumValue_ = maximumIndependentVariable1_;

        break;

    case 2:

        maximumValue_ = maximumIndependentVariable2_;

        break;

    default:

        maximumValue_ = -0.0;

        cerr << " Only 2 independent variables, variable "
             <<  parameterIndex << " does not exist when"
             << "getting maximum value, returning -0.0" << endl;
    }

    // Return maximum value.
    return maximumValue_;
}

//! Apply transformation to vehicle part.
void SingleSurfaceGeometry::transformPoint( VectorXd& point )
{
    // Apply scaling, rotation and translation operations.
    point = scalingMatrix_ * point;
    point = rotationMatrix_ * point;
    point = point + offset_;
}

}

// End of file.

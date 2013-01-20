/*    Copyright (c) 2010-2013, Delft University of Technology
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
 *      120323    D. Dirkx          Removed set functions; moved functionality to constructor.
 *
 *    References
 *
 *    Notes
 *      Contents of this file used to be in singleGeometry.cpp, but as this class has been split
 *      into single and composite surface geometry, the contents have been moved, with most of the
 *      SurfaceGeometry class now belonging to the SingleSurfaceGeometry class.
 *
 */

#include <iostream>

#include "Tudat/Mathematics/GeometricShapes/singleSurfaceGeometry.h"

namespace tudat
{
namespace geometric_shapes
{

using std::cerr;
using std::endl;

//! Set minimum value of independent variable.
void SingleSurfaceGeometry::setMinimumIndependentVariable( const int parameterIndex,
                                                           const double minimumValue )
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
        const int parameterIndex, const double maximumValue )
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
double SingleSurfaceGeometry::getMinimumIndependentVariable( const int parameterIndex )
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
double SingleSurfaceGeometry::getMaximumIndependentVariable( const int parameterIndex )
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

        cerr << "Only 2 independent variables, variable " << parameterIndex << " does not exist "
             << "when getting maximum value, returning -0.0" << endl;
    }

    // Return maximum value.
    return maximumValue_;
}

//! Apply transformation to vehicle part.
void SingleSurfaceGeometry::transformPoint( Eigen::VectorXd& point )
{
    // Apply scaling, rotation and translation operations.
    point = scalingMatrix_ * point;
    point = rotationMatrix_ * point;
    point = point + offset_;
}

} // namespace geometric_shapes
} // namespace tudat

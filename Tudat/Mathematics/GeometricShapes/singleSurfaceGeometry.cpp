/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      Contents of this file used to be in singleGeometry.cpp, but as this class has been split
 *      into single and composite surface geometry, the contents have been moved, with most of the
 *      SurfaceGeometry class now belonging to the SingleSurfaceGeometry class.
 *
 */

#include "Tudat/Mathematics/GeometricShapes/singleSurfaceGeometry.h"

namespace tudat
{
namespace geometric_shapes
{

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

        std::string errorMessage =  "Only 2 independent variables, variable " +
                std::to_string( parameterIndex ) + " does not exist when setting minimum value";
        throw std::runtime_error( errorMessage );
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

        std::string errorMessage =  "Only 2 independent variables, variable " +
                std::to_string( parameterIndex ) + " does not exist when setting maximum value";
        throw std::runtime_error( errorMessage );
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

        std::string errorMessage =  "Only 2 independent variables, variable " +
                std::to_string( parameterIndex ) + " does not exist when getting minimum value";
        throw std::runtime_error( errorMessage );
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
        std::string errorMessage =  "Only 2 independent variables, variable " +
                std::to_string( parameterIndex ) + " does not exist when getting maximum value";
        throw std::runtime_error( errorMessage );
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

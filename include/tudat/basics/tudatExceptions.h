/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_EXCEPTIONS_H
#define TUDAT_EXCEPTIONS_H

#include <iostream>
#include <string>
#include <boost/lexical_cast.hpp>

namespace tudat
{

namespace exceptions
{

/// @brief Generic exception class for Tudat.
///
/// This class is a generic base class for all Tudat exceptions.
/// It inherits from std::runtime_error, which makes it backwards compatible with
/// the previous behavior, where only std::runtime_error was used.
class TudatError : public std::runtime_error
{
private:
public:
    /// @brief Default constructor for TudatError.
    /// @param errorMessage Error message to be displayed.
    TudatError( const std::string& errorMessage ): std::runtime_error( errorMessage ) { }
    ~TudatError( ) { }
};

/// @brief Error class when an interpolator is requested to return a value outside of its bounds.
/// @tparam T Type of the independent variable data.
template< typename T >
class InterpolationOutOfBoundsError : public TudatError
{
private:
public:
    /// @brief Constructor for one-dimensional interpolation out of bounds error.
    /// @param requestedValue Requested value which is out of bounds.
    /// @param lowerBound Value of the lower bound of the interpolation data.
    /// @param upperBound Value of the upper bound of the interpolation data.
    InterpolationOutOfBoundsError( const T requestedValue, const T lowerBound, const T upperBound ):
        exceptions::TudatError( "Error in interpolator, requesting data point outside of boundaries, requested data at " +
                                boost::lexical_cast< std::string >( requestedValue ) + " but limit values are " +
                                boost::lexical_cast< std::string >( lowerBound ) + " and " +
                                boost::lexical_cast< std::string >( upperBound ) ),
        requestedValue( requestedValue ), lowerBound( lowerBound ), upperBound( upperBound )
    { }
    /// @brief Constructor for multi-dimensional interpolation out of bounds error.
    /// @param requestedValue Requested value which is out of bounds.
    /// @param lowerBound Value of the lower bound of the interpolation data in the requested dimension.
    /// @param upperBound Value of the upper bound of the interpolation data in the requested dimension.
    /// @param requestDimension Dimension in which the requested data is out of bounds.
    InterpolationOutOfBoundsError( const T requestedValue, const T lowerBound, const T upperBound, const unsigned int requestDimension ):
        exceptions::TudatError( "Error in interpolator, requesting data point outside of boundaries, requested data of dimension " +
                                boost::lexical_cast< std::string >( requestDimension ) + " at " +
                                boost::lexical_cast< std::string >( requestedValue ) + " but limit values are " +
                                boost::lexical_cast< std::string >( lowerBound ) + " and " +
                                boost::lexical_cast< std::string >( upperBound ) ),
        requestedValue( requestedValue ), lowerBound( lowerBound ), upperBound( upperBound )
    { }

    ~InterpolationOutOfBoundsError( ) { }

    const T requestedValue;
    const T lowerBound;
    const T upperBound;
};

/// @brief Error class when a Lagrange interpolator is requested to return a value outside of its reliable bounds.
/// @tparam T Type of the independent variable data.
template< typename T >
class LagrangeInterpolationOutOfBoundsError : public InterpolationOutOfBoundsError< T >
{
private:
public:
    /// @brief Constructor for Lagrange interpolation out of bounds error.
    /// @param requestedValue Requested value which is out of the reliable bounds of the Lagrange interpolation.
    /// @param lowerBound Value of the lower reliable bound of the interpolation data.
    /// @param upperBound Value of the upper reliable bound of the interpolation data.
    LagrangeInterpolationOutOfBoundsError( const T requestedValue, const T lowerBound, const T upperBound ):
        exceptions::InterpolationOutOfBoundsError< T >( requestedValue, lowerBound, upperBound )
    { }

    ~LagrangeInterpolationOutOfBoundsError( ) { }
};

/// @brief Error class for the number of maximum iterations have been exceeded.
class MaximumIterationsExceededError : public TudatError
{
private:
public:
    /// @brief Constructor for maximum iterations exceeded error.
    /// @param numberOfIterations Number of iterations that have been performed, which exceeded the maximum.
    /// @param maximumNumberOfIterations Maximum number of iterations that were allowed.
    MaximumIterationsExceededError( const unsigned int numberOfIterations, const unsigned int maximumNumberOfIterations ):
        exceptions::TudatError( "Root-finder did not converge within maximum number of iterations! " +
                                std::to_string( numberOfIterations ) + " iterations completed, maximum number of iterations is " +
                                std::to_string( maximumNumberOfIterations ) ),
        numberOfIterations( numberOfIterations ), maximumNumberOfIterations( maximumNumberOfIterations )
    { }
    ~MaximumIterationsExceededError( ) { }

    const unsigned int numberOfIterations;
    const unsigned int maximumNumberOfIterations;
};

/// @brief Generic error class when a step size violation occurs that cannot be resolved.
class StepSizeViolationError : public TudatError
{
private:
public:
    StepSizeViolationError( const std::string& errorMessage ): exceptions::TudatError( errorMessage ) { }
    ~StepSizeViolationError( ) { }
};

/// @brief Error class when the requested step size is smaller than the minimum step size.
/// @tparam TimeStepType Type of the time step.
template< typename TimeStepType >
class MinimumStepSizeViolatedError : public StepSizeViolationError
{
private:
public:
    /// @brief Constructor for minimum step size violated error.
    /// @param minimumStepSize Minimum step size that is violated.
    /// @param recommendedStepSize Step size recommended by the step-size controller which violates the minimumStepSize.
    MinimumStepSizeViolatedError( const TimeStepType minimumStepSize, const TimeStepType recommendedStepSize ):
        exceptions::StepSizeViolationError( "Error in step-size control, minimum step size " + std::to_string( minimumStepSize ) +
                                            " is higher than required time step " + std::to_string( recommendedStepSize ) ),
        minimumStepSize( minimumStepSize ), recommendedStepSize( recommendedStepSize )
    { }
    ~MinimumStepSizeViolatedError( ) { }

    const TimeStepType minimumStepSize;
    const TimeStepType recommendedStepSize;
};

}  // namespace exceptions

}  // namespace tudat

#endif  // TUDAT_EXCEPTIONS_H

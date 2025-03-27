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

class TudatError : public std::runtime_error
{
private:
public:
    TudatError( const std::string& errorMessage ): std::runtime_error( errorMessage ) { }
    ~TudatError( ) { }
};

template< typename T >
class InterpolationOutOfBoundsError : public TudatError
{
private:
public:
    InterpolationOutOfBoundsError( const T requestedValue, const T lowerBound, const T upperBound ):
        exceptions::TudatError( "Error in interpolator, requesting data point outside of boundaries, requested data at " +
                                boost::lexical_cast< std::string >( requestedValue ) + " but limit values are " +
                                boost::lexical_cast< std::string >( lowerBound ) + " and " +
                                boost::lexical_cast< std::string >( upperBound ) ),
        requestedValue( requestedValue ), lowerBound( lowerBound ), upperBound( upperBound )
    { }
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

template< typename T >
class LagrangeInterpolationOutOfBoundsError : public InterpolationOutOfBoundsError< T >
{
private:
public:
    LagrangeInterpolationOutOfBoundsError( const T requestedValue, const T lowerBound, const T upperBound ):
        exceptions::InterpolationOutOfBoundsError< T >( requestedValue, lowerBound, upperBound )
    { }

    ~LagrangeInterpolationOutOfBoundsError( ) { }
};

class MaximumIterationsExceededError : public TudatError
{
private:
public:
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

class StepSizeViolationError : public TudatError
{
private:
public:
    StepSizeViolationError( const std::string& errorMessage ): exceptions::TudatError( errorMessage ) { }
    ~StepSizeViolationError( ) { }
};

template< typename TimeStepType >
class MinimumStepSizeViolatedError : public StepSizeViolationError
{
private:
public:
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

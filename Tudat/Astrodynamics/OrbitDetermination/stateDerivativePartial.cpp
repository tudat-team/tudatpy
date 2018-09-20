/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/stateDerivativePartial.h"

namespace tudat
{

namespace orbit_determination
{

//! Function to evaluate the negative value of a parameter partial.
void evaluateNegativeParameterPartialFunction(
        const std::function< void( Eigen::MatrixXd& ) > parameterPartialFunction,
        Eigen::MatrixXd& partial )
{
    parameterPartialFunction( partial );
    partial *= -1.0;
}

//! Function to evaluate the subtraction of two parameter partials.
void evaluateSubtractedParameterPartialFunction(
        const std::function< void( Eigen::MatrixXd& ) > firstParameterPartialFunction,
        const std::function< void( Eigen::MatrixXd& ) > parameterPartialFunctionToSubtract,
        Eigen::MatrixXd& partial )
{
    firstParameterPartialFunction( partial );

    Eigen::MatrixXd subtractedPartial = Eigen::MatrixXd::Zero( partial.rows( ), partial.cols( ) );
    parameterPartialFunctionToSubtract( subtractedPartial );

    partial -= subtractedPartial;
}

//! Function to evaluate the addition of two parameter partials.
void evaluateAddedParameterPartialFunction(
        const std::function< void( Eigen::MatrixXd& ) > firstParameterPartialFunction,
        const std::function< void( Eigen::MatrixXd& ) > parameterPartialFunctionToAdd,
        Eigen::MatrixXd& partial )
{
    firstParameterPartialFunction( partial );

    Eigen::MatrixXd addedPartial = Eigen::MatrixXd::Zero( partial.rows( ), partial.cols( ) );
    parameterPartialFunctionToAdd( addedPartial );

    partial += addedPartial;
}

//! Create a parameter partial function obtained from the subtraction of two such function results.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > createMergedParameterPartialFunction(
        const std::pair< std::function< void( Eigen::MatrixXd& ) >, int >& partialFunctionOfAccelerationToAdd,
        const std::pair< std::function< void( Eigen::MatrixXd& ) >, int >& partialFunctionOfAccelerationToSubtract )
{
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int >  parameterPartialFunction;

    // Check partial size and act accordingly.
    if( ( partialFunctionOfAccelerationToAdd.second == 0 ) &&
            ( partialFunctionOfAccelerationToSubtract.second == 0 ) )
    {
        parameterPartialFunction = std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );
    }
    else if( partialFunctionOfAccelerationToSubtract.second == 0 )
    {
        parameterPartialFunction = partialFunctionOfAccelerationToAdd;
    }
    else if( partialFunctionOfAccelerationToAdd.second == 0 )
    {
        parameterPartialFunction = std::make_pair( std::bind(
                                                       &evaluateNegativeParameterPartialFunction,
                                                       partialFunctionOfAccelerationToSubtract.first, std::placeholders::_1 ),
                                                   partialFunctionOfAccelerationToSubtract.second );
    }
    // Partial size must be equal if both non-zero
    else if( partialFunctionOfAccelerationToSubtract.second !=
             partialFunctionOfAccelerationToAdd.second )
    {
        throw std::runtime_error( "Error when making merged parameter partial function, separate functions are incompatible" );
    }
    else
    {
        parameterPartialFunction = std::make_pair( std::bind( &evaluateSubtractedParameterPartialFunction,
                                                               partialFunctionOfAccelerationToAdd.first,
                                                               partialFunctionOfAccelerationToSubtract.first, std::placeholders::_1 ),
                                                  partialFunctionOfAccelerationToSubtract.second );
    }
    return parameterPartialFunction;
}

//! Function to create a parameter partial evaluation function, obtained by adding or subtracting a given partial
//! w.r.t. a double parameter from 2 state derivative partial models.
std::function< void( Eigen::MatrixXd& ) > getCombinedCurrentDoubleParameterFunction(
        const std::shared_ptr< StateDerivativePartial > firstPartial,
        const std::shared_ptr< StateDerivativePartial > secondPartial,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameterObject,
        const int firstPartialSize, const int secondPartialSize,
        const bool subtractPartials )
{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;

    // Get two partial functions.
    std::function< void( Eigen::MatrixXd& ) > firstPartialFunction =
            std::bind( &StateDerivativePartial::getCurrentDoubleParameterPartial, firstPartial, parameterObject, std::placeholders::_1 );
    std::function< void( Eigen::MatrixXd& ) > secondPartialFunction =
            std::bind( &StateDerivativePartial::getCurrentDoubleParameterPartial, secondPartial, parameterObject, std::placeholders::_1 );

    // If both partial function sizes are zero, cannot create partial.
    if( firstPartialSize == 0 && secondPartialSize == 0 )
    {
        throw std::runtime_error( "Error when getting combined current partial size, both partials have size zero " );
    }
    // Check other cases and combine partial functions accordingly
    else if( secondPartialSize == 0 )
    {
        partialFunction = firstPartialFunction;
    }
    else if( firstPartialSize == 0 )
    {
        if( subtractPartials )
        {
            partialFunction = std::bind( &evaluateNegativeParameterPartialFunction,
                                           secondPartialFunction, std::placeholders::_1 );
        }
        else
        {
            partialFunction = secondPartialFunction;
        }
    }
    else if( firstPartialSize == secondPartialSize )
    {
        if( subtractPartials )
        {
            partialFunction = std::bind( &evaluateSubtractedParameterPartialFunction,
                                           firstPartialFunction,
                                           secondPartialFunction, std::placeholders::_1 );
        }
        else
        {
            partialFunction = std::bind( &evaluateAddedParameterPartialFunction,
                                           firstPartialFunction,
                                           secondPartialFunction, std::placeholders::_1 );
        }
    }
    // Partial size must be equal if both non-zero
    else
    {
        throw std::runtime_error(
                    "Error when getting combined current partial size, both partials have different non-zero size." );
    }
    return partialFunction;
}

//! Function to create a parameter partial evaluation function, obtained by adding or subtracting a given partial
//! w.r.t. a vector parameter from 2 state derivative partial models.
std::function< void( Eigen::MatrixXd& ) > getCombinedCurrentVectorParameterFunction(
        const std::shared_ptr< StateDerivativePartial > firstPartial,
        const std::shared_ptr< StateDerivativePartial > secondPartial,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameterObject,
        const int firstPartialSize, const int secondPartialSize,
        const bool subtractPartials )
{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;

    // Get two partial functions.
    std::function< void( Eigen::MatrixXd& ) > firstPartialFunction =
            std::bind( &StateDerivativePartial::getCurrentVectorParameterPartial, firstPartial, parameterObject, std::placeholders::_1 );
    std::function< void( Eigen::MatrixXd& ) > secondPartialFunction =
            std::bind( &StateDerivativePartial::getCurrentVectorParameterPartial, secondPartial, parameterObject, std::placeholders::_1 );

    // If both partial function sizes are zero, cannot create partial.
    if( firstPartialSize == 0 && secondPartialSize == 0 )
    {
        throw std::runtime_error( "Error when getting combined current partial size, both partials have size zero " );
    }
    // Check other cases and combine partial functions accordingly
    else if( secondPartialSize == 0 )
    {
        partialFunction = firstPartialFunction;
    }
    else if( firstPartialSize == 0 )
    {
        if( subtractPartials )
        {
            partialFunction = std::bind( &evaluateNegativeParameterPartialFunction,
                                           secondPartialFunction, std::placeholders::_1 );
        }
        else
        {
            partialFunction = secondPartialFunction;

        }
    }
    else if( firstPartialSize == secondPartialSize )
    {
        if( subtractPartials )
        {
            partialFunction = std::bind( &evaluateSubtractedParameterPartialFunction,
                                           firstPartialFunction,
                                           secondPartialFunction, std::placeholders::_1 );
        }
        else
        {
            partialFunction = std::bind( &evaluateAddedParameterPartialFunction,
                                           firstPartialFunction,
                                           secondPartialFunction, std::placeholders::_1 );
        }

    }
    // Partial size must be equal if both non-zero
    else
    {
        throw std::runtime_error( "Error when getting combined current partial size, both partials have different non-zero size." );
    }
    return partialFunction;
}

}

}


/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/mutualSphericalHarmonicGravityPartial.h"

namespace tudat
{

namespace acceleration_partials
{

//! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > MutualSphericalHarmonicsGravityPartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > parameterPartial;

    // Check if parameter is gravitational parameter of either
    if( parameter->getParameterName( ).first == estimatable_parameters::gravitational_parameter )
    {
        // Check of acceleration depends on gravitational parameter of given body
        if( parameter->getParameterName( ).second.first == acceleratingBody_ ||
                ( parameter->getParameterName( ).second.first == acceleratedBody_ && accelerationUsesMutualAttraction_ ) )
        {
            parameterPartial = std::make_pair(
                        std::bind( &MutualSphericalHarmonicsGravityPartial::wrtGravitationalParameter, this, std::placeholders::_1 ), 1 );
        }
        else
        {
            parameterPartial = std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );
        }
    }
    else
    {
        // Get partial functions for constituent partial objects
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromBodyExertingAcceleration =
                accelerationPartialOfShExpansionOfBodyExertingAcceleration_->getParameterPartialFunction( parameter );
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromBodyUndergoingAcceleration =
                accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->getParameterPartialFunction( parameter );

        // Combine partial functions
        parameterPartial = orbit_determination::createMergedParameterPartialFunction(
                    partialFunctionFromBodyExertingAcceleration, partialFunctionFromBodyUndergoingAcceleration );
    }
    return parameterPartial;
}

//! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > MutualSphericalHarmonicsGravityPartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromBodyExertingAcceleration =
            accelerationPartialOfShExpansionOfBodyExertingAcceleration_->getParameterPartialFunction( parameter );
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromBodyUndergoingAcceleration =
            accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->getParameterPartialFunction( parameter );

    return  orbit_determination::createMergedParameterPartialFunction(
                partialFunctionFromBodyExertingAcceleration, partialFunctionFromBodyUndergoingAcceleration );
}

//! Function to set a dependency of this partial object w.r.t. a given double parameter.
int MutualSphericalHarmonicsGravityPartial::setParameterPartialUpdateFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    int partialSize = 0;

    // Check if parameter is gravitational parameter of either
    if( parameter->getParameterName( ).first == estimatable_parameters::gravitational_parameter )
    {
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int >  partialFunction  =
                getParameterPartialFunction( parameter );
        partialSize = partialFunction.second;

        if( partialFunction.second > 0 )
        {
            parameterDoublePartialFunctions_[ parameter ] = partialFunction.first;
            isCurrentDoubleParameterPartialSet_[ parameter ] = 0;
            currentDoubleParameterPartials_[ parameter ] = Eigen::MatrixXd( accelerationSize_, 1 );
        }
    }
    else
    {
        // Get partial function for acceleration from body exerting acceleration
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromExertingExpansion =
                accelerationPartialOfShExpansionOfBodyExertingAcceleration_->getParameterPartialFunction( parameter );
        if( partialFunctionFromExertingExpansion.second > 0 )
        {
            accelerationPartialOfShExpansionOfBodyExertingAcceleration_->setParameterPartialUpdateFunction( parameter );
        }

        // Get partial function for acceleration from body undergoing acceleration
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromUndergoingGravity =
                accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->getParameterPartialFunction( parameter );
        if( partialFunctionFromUndergoingGravity.second > 0 )
        {
            accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->setParameterPartialUpdateFunction( parameter );
        }

        // Update this partial object with new dependencies
        if( partialFunctionFromUndergoingGravity.second > 0 || partialFunctionFromExertingExpansion.second > 0 )
        {
            parameterDoublePartialFunctions_[ parameter ] =
                    getCombinedCurrentDoubleParameterFunction(
                        accelerationPartialOfShExpansionOfBodyExertingAcceleration_,
                        accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_,
                        parameter, partialFunctionFromExertingExpansion.second, partialFunctionFromUndergoingGravity.second, 1 );
            isCurrentDoubleParameterPartialSet_[ parameter ] = 0;
            currentDoubleParameterPartials_[ parameter ] = Eigen::MatrixXd( accelerationSize_, 1 );
        }
        partialSize = std::max( partialFunctionFromExertingExpansion.second, partialFunctionFromUndergoingGravity.second );
    }
    return partialSize;
}

//! Function to set a dependency of this partial object w.r.t. a given vector parameter.
int MutualSphericalHarmonicsGravityPartial::setParameterPartialUpdateFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    // Get partial function for acceleration from body exerting acceleration
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromExertingExpansion =
            accelerationPartialOfShExpansionOfBodyExertingAcceleration_->getParameterPartialFunction( parameter );
    if( partialFunctionFromExertingExpansion.second > 0 )
    {
        accelerationPartialOfShExpansionOfBodyExertingAcceleration_->setParameterPartialUpdateFunction( parameter );
    }

    // Get partial function for acceleration from body undergoing acceleration
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromUndergoingGravity =
            accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->getParameterPartialFunction( parameter );
    if( partialFunctionFromUndergoingGravity.second > 0 )
    {
        accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->setParameterPartialUpdateFunction( parameter );
    }

    // Update this partial object with new dependencies
    if( partialFunctionFromUndergoingGravity.second > 0 || partialFunctionFromExertingExpansion.second > 0 )
    {
        parameterVectorPartialFunctions_[ parameter ] =
                getCombinedCurrentVectorParameterFunction(
                    accelerationPartialOfShExpansionOfBodyExertingAcceleration_,
                    accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_,
                    parameter, partialFunctionFromExertingExpansion.second, partialFunctionFromUndergoingGravity.second, 1 );
        isCurrentVectorParameterPartialSet_[ parameter ] = 0;
        currentVectorParameterPartials_[ parameter ] = Eigen::MatrixXd( accelerationSize_, parameter->getParameterSize( ) );

    }
    return std::max( partialFunctionFromExertingExpansion.second, partialFunctionFromUndergoingGravity.second );
}

}

}

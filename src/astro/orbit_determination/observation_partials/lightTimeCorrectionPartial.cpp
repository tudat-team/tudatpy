/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */





#include "tudat/astro/orbit_determination/observation_partials/lightTimeCorrectionPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/firstOrderRelativisticPartial.h"

namespace tudat
{

namespace observation_partials
{

//! Function to get the function returning the light-time correction partial for given correction partial and parameter.

std::pair< std::function< LightTimeCorrectionPartial::SingleOneWayRangePartialReturnType(
        const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) >, bool >
getLightTimeParameterPartialFunction(
        const estimatable_parameters::EstimatebleParameterIdentifier parameterId,
        const std::shared_ptr< LightTimeCorrectionPartial > lightTimeCorrectionPartial )
{
    // Declare return type, set second part to 0 (no dependency found).
    std::pair< std::function< LightTimeCorrectionPartial::SingleOneWayRangePartialReturnType(
                const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) >, bool > partialFunction;
    partialFunction.second = 0;

    // Check type of light-time correction
    switch( lightTimeCorrectionPartial->getCorrectionType( ) )
    {
    // Correction type of 1st-order relativistic
    case observation_models::first_order_relativistic:
    {
        // Check consistency of input.
        std::shared_ptr< FirstOrderRelativisticLightTimeCorrectionPartial > currentLightTimeCorrectorPartial =
                std::dynamic_pointer_cast< FirstOrderRelativisticLightTimeCorrectionPartial >( lightTimeCorrectionPartial );
        if( currentLightTimeCorrectorPartial == nullptr )
        {
            std::string errorMessage = "Error when getting light time correction partial function, type " +
                    std::to_string( lightTimeCorrectionPartial->getCorrectionType( ) ) +
                       "is inconsistent.";
            throw std::runtime_error( errorMessage );
        }

        // Set partial of gravitational parameter
        if( parameterId.first == estimatable_parameters::gravitational_parameter )
        {
            // Retrieve function from FirstOrderRelativisticLightTimeCorrectionPartial if correction depends on
            // body associated with parameter.
            std::vector< std::string > perturbingBodies = currentLightTimeCorrectorPartial->getPerturbingBodies( );
            std::vector< std::string >::iterator findIterator = std::find(
                        perturbingBodies.begin( ), perturbingBodies.end( ), parameterId.second.first );
            if( findIterator != perturbingBodies.end( ) )
            {
                int bodyIndex = std::distance( perturbingBodies.begin( ),  findIterator );
                partialFunction = std::make_pair(
                            std::bind( &FirstOrderRelativisticLightTimeCorrectionPartial::wrtBodyGravitationalParameter,
                                         currentLightTimeCorrectorPartial, std::placeholders::_1, std::placeholders::_2, bodyIndex ), 1 );
            }
        }
        else if( parameterId.first == estimatable_parameters::ppn_parameter_gamma )
        {
            partialFunction = std::make_pair(
                        std::bind( &FirstOrderRelativisticLightTimeCorrectionPartial::wrtPpnParameterGamma,
                                     currentLightTimeCorrectorPartial, std::placeholders::_1, std::placeholders::_2 ), 1 );
        }
        break;
    }
    default:
        std::string errorMessage = "Error, light time correction type " + std::to_string(
                    lightTimeCorrectionPartial->getCorrectionType( ) ) + "not found when creating partial ";
        throw std::runtime_error( errorMessage );
    }

    return partialFunction;
}

std::pair< std::function< LightTimeCorrectionPartial::SingleOneWayRangeGradientPartialReturnType(
    const std::vector< Eigen::Vector6d >&, const std::vector< double >&, const observation_models::LinkEndType ) >, bool >
getLightTimeGradientParameterPartialFunction(
    const estimatable_parameters::EstimatebleParameterIdentifier parameterId,
    const std::shared_ptr< LightTimeCorrectionPartial > lightTimeCorrectionPartial )
{
    // Declare return type, set second part to 0 (no dependency found).
    std::pair< std::function< LightTimeCorrectionPartial::SingleOneWayRangeGradientPartialReturnType(
        const std::vector< Eigen::Vector6d >&, const std::vector< double >&, const observation_models::LinkEndType ) >, bool > partialFunction;
    partialFunction.second = 0;

    // Check type of light-time correction
    switch( lightTimeCorrectionPartial->getCorrectionType( ) )
    {
    // Correction type of 1st-order relativistic
    case observation_models::first_order_relativistic:
    {
        // Check consistency of input.
        std::shared_ptr< FirstOrderRelativisticLightTimeCorrectionPartial > currentLightTimeCorrectorPartial =
            std::dynamic_pointer_cast< FirstOrderRelativisticLightTimeCorrectionPartial >( lightTimeCorrectionPartial );
        if( currentLightTimeCorrectorPartial == nullptr )
        {
            std::string errorMessage = "Error when getting light time correction partial function, type " +
                                       std::to_string( lightTimeCorrectionPartial->getCorrectionType( ) ) +
                                       "is inconsistent.";
            throw std::runtime_error( errorMessage );
        }

        // Set partial of gravitational parameter
        if( parameterId.first == estimatable_parameters::ppn_parameter_gamma )
        {
            partialFunction = std::make_pair(
                std::bind( &FirstOrderRelativisticLightTimeCorrectionPartial::gradientWrtPpnParameterGamma,
                           currentLightTimeCorrectorPartial, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3 ), 1 );
        }
        break;
    }
    default:
        std::string errorMessage = "Error, light time correction type " + std::to_string(
            lightTimeCorrectionPartial->getCorrectionType( ) ) + "not found when creating partial ";
        throw std::runtime_error( errorMessage );
    }

    return partialFunction;
}

}

}

/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#ifndef TUDAT_CREATEOBSERVATIONMANAGER_H
#define TUDAT_CREATEOBSERVATIONMANAGER_H

#include <vector>
#include <map>

#include <memory>

#include "tudat/astro/observation_models/observationManager.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/observationBiasParameter.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/estimation_setup/createObservationPartials.h"
#include "tudat/astro/observation_models/oneWayRangeObservationModel.h"
#include "tudat/simulation/propagation_setup/dependentVariablesInterface.h"

namespace tudat
{

namespace observation_models
{

//! Function to perform the closure a single observation bias and a single estimated bias parameter.
/*!
 *  Function to perform the closure a single observation bias and a single estimated bias parameter. Estimated parameter objects
 *  are typically created prior to observation models. This function must be called for the estimated parameter object creation
 *  to be finalized, in the case link properties are estimated (e.g. observation biases).
 */
template< int ObservationSize = 1 >
void performObservationParameterEstimationClosureForSingleModelSet(
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
        const std::shared_ptr< ObservationBias< ObservationSize > > observationBias,
        const LinkEnds linkEnds,
        const ObservableType observableType )
{
    ObservationBiasTypes biasType = getObservationBiasType( observationBias );

    // Check if bias type is multi-type
    if( biasType == multiple_observation_biases )
    {
        // Test input consistency
        std::shared_ptr< MultiTypeObservationBias< ObservationSize > > multiTypeBias =
                std::dynamic_pointer_cast< MultiTypeObservationBias< ObservationSize > >( observationBias );
        if( multiTypeBias == nullptr )
        {
            throw std::runtime_error( "Error, cannot perform bias closure for multi-type bias, inconsistent bias types" );
        }

        // Perform closure for each constituent bias object
        for( unsigned int i = 0; i < multiTypeBias->getBiasList( ).size( ); i++ )
        {
            performObservationParameterEstimationClosureForSingleModelSet(
                        parameter, multiTypeBias->getBiasList( ).at( i ), linkEnds, observableType );
        }
    }
    else
    {
        // Check bias type
        switch( parameter->getParameterName( ).first )
        {
        case estimatable_parameters::constant_additive_observation_bias:
        {
            // Test input consistency
            std::shared_ptr< estimatable_parameters::ConstantObservationBiasParameter > biasParameter =
                    std::dynamic_pointer_cast< estimatable_parameters::ConstantObservationBiasParameter >(
                        parameter );
            if( biasParameter == nullptr )
            {
                throw std::runtime_error( "Error, cannot perform bias closure for additive bias, inconsistent bias types" );
            }

            // Check if bias object is of same type as estimated parameter
            std::shared_ptr< ConstantObservationBias< ObservationSize > > constantBiasObject =
                    std::dynamic_pointer_cast< ConstantObservationBias< ObservationSize > >( observationBias );
            if( constantBiasObject != nullptr )
            {
                // Check if bias and parameter link properties are equal
                if( linkEnds == biasParameter->getLinkEnds( ) &&
                        observableType == biasParameter->getObservableType( ) )
                {
                    biasParameter->setObservationBiasFunctions(
                                std::bind( &ConstantObservationBias< ObservationSize >::getTemplateFreeConstantObservationBias,
                                           constantBiasObject ),
                                std::bind( &ConstantObservationBias< ObservationSize >::resetConstantObservationBiasTemplateFree,
                                           constantBiasObject, std::placeholders::_1 ) );
                }
            }
            break;
        }
        case estimatable_parameters::arcwise_constant_additive_observation_bias:
        {

            // Test input consistency
            std::shared_ptr< estimatable_parameters::ArcWiseObservationBiasParameter > biasParameter =
                    std::dynamic_pointer_cast< estimatable_parameters::ArcWiseObservationBiasParameter >(
                        parameter );
            if( biasParameter == nullptr )
            {
                throw std::runtime_error( "Error, cannot perform bias closure for arc-wise additive bias, inconsistent bias types" );
            }

            // Check if bias object is of same type as estimated parameter
            std::shared_ptr< ConstantArcWiseObservationBias< ObservationSize > > constantBiasObject =
                    std::dynamic_pointer_cast< ConstantArcWiseObservationBias< ObservationSize > >( observationBias );
            if( constantBiasObject != nullptr )
            {
                // Check if bias and parameter link properties are equal
                if( ( linkEnds == biasParameter->getLinkEnds( ) ) &&
                        ( observableType == biasParameter->getObservableType( ) ) &&
                        ( constantBiasObject->getLinkEndIndexForTime( ) == biasParameter->getLinkEndIndex( ) ) &&
                        ( constantBiasObject->getArcStartTimes( ).size( ) == biasParameter->getArcStartTimes( ).size( ) ) )
                {
                    bool doTimesMatch = true;
                    for( unsigned int i = 0; i < constantBiasObject->getArcStartTimes( ).size( ); i++ )
                    {
                        if( std::fabs( constantBiasObject->getArcStartTimes( ).at( i ) -
                                       biasParameter->getArcStartTimes( ).at( i ) ) >  std::max(
                                    1.0E-15 * std::fabs( constantBiasObject->getArcStartTimes( ).at( i ) ),
                                    1.0E-15 * std::fabs( biasParameter->getArcStartTimes( ).at( i ) ) ) )
                        {
                            doTimesMatch = false;
                        }
                    }

                    if( doTimesMatch == true )
                    {
                        biasParameter->setObservationBiasFunctions(
                                    std::bind( &ConstantArcWiseObservationBias< ObservationSize >::getTemplateFreeConstantObservationBias,
                                               constantBiasObject ),
                                    std::bind( &ConstantArcWiseObservationBias< ObservationSize >::resetConstantObservationBiasTemplateFree,
                                               constantBiasObject, std::placeholders::_1 ) );
                        biasParameter->setLookupScheme( constantBiasObject->getLookupScheme( ) );
                    }
                }
            }
            break;
        }
        case estimatable_parameters::constant_relative_observation_bias:
        {
            // Test input consistency
            std::shared_ptr< estimatable_parameters::ConstantObservationBiasParameter > biasParameter =
                    std::dynamic_pointer_cast< estimatable_parameters::ConstantObservationBiasParameter >(
                        parameter );
            if( biasParameter == nullptr )
            {
                throw std::runtime_error( "Error, cannot perform bias closure for additive bias, inconsistent bias types" );
            }

            // Check if bias object is of same type as estimated parameter
            std::shared_ptr< ConstantRelativeObservationBias< ObservationSize > > constantBiasObject =
                    std::dynamic_pointer_cast< ConstantRelativeObservationBias< ObservationSize > >( observationBias );
            if( constantBiasObject != nullptr )
            {
                // Check if bias and parameter link properties are equal
                if( linkEnds == biasParameter->getLinkEnds( ) &&
                        observableType == biasParameter->getObservableType( ) )
                {
                    biasParameter->setObservationBiasFunctions(
                                std::bind( &ConstantRelativeObservationBias< ObservationSize >::getTemplateFreeConstantObservationBias,
                                           constantBiasObject ),
                                std::bind( &ConstantRelativeObservationBias< ObservationSize >::resetConstantObservationBiasTemplateFree,
                                           constantBiasObject, std::placeholders::_1 ) );
                }
            }
            break;
        }
        case estimatable_parameters::arcwise_constant_relative_observation_bias:
        {

            // Test input consistency
            std::shared_ptr< estimatable_parameters::ArcWiseObservationBiasParameter > biasParameter =
                    std::dynamic_pointer_cast< estimatable_parameters::ArcWiseObservationBiasParameter >(
                        parameter );
            if( biasParameter == nullptr )
            {
                throw std::runtime_error( "Error, cannot perform bias closure for arc-wise relative bias, inconsistent bias types" );
            }

            // Check if bias object is of same type as estimated parameter
            std::shared_ptr< ConstantRelativeArcWiseObservationBias< ObservationSize > > constantBiasObject =
                    std::dynamic_pointer_cast< ConstantRelativeArcWiseObservationBias< ObservationSize > >( observationBias );
            if( constantBiasObject != nullptr )
            {
                // Check if bias and parameter link properties are equal
                if( ( linkEnds == biasParameter->getLinkEnds( ) ) &&
                        ( observableType == biasParameter->getObservableType( ) ) &&
                        ( constantBiasObject->getLinkEndIndexForTime( ) == biasParameter->getLinkEndIndex( ) ) &&
                        ( constantBiasObject->getArcStartTimes( ).size( ) == biasParameter->getArcStartTimes( ).size( ) ) )
                {
                    bool doTimesMatch = true;
                    for( unsigned int i = 0; i < constantBiasObject->getArcStartTimes( ).size( ); i++ )
                    {
                        if( std::fabs( constantBiasObject->getArcStartTimes( ).at( i ) -
                                       biasParameter->getArcStartTimes( ).at( i ) ) >  std::max(
                                    1.0E-15 * std::fabs( constantBiasObject->getArcStartTimes( ).at( i ) ),
                                    1.0E-15 * std::fabs( biasParameter->getArcStartTimes( ).at( i ) ) ) )
                        {
                            doTimesMatch = false;
                        }
                    }

                    if( doTimesMatch == true )
                    {
                        biasParameter->setObservationBiasFunctions(
                                    std::bind( &ConstantRelativeArcWiseObservationBias< ObservationSize >::getTemplateFreeConstantObservationBias,
                                               constantBiasObject ),
                                    std::bind( &ConstantRelativeArcWiseObservationBias< ObservationSize >::resetConstantObservationBiasTemplateFree,
                                               constantBiasObject, std::placeholders::_1 ) );
                        biasParameter->setLookupScheme( constantBiasObject->getLookupScheme( ) );
                    }
                }
            }
            break;
        }
        case estimatable_parameters::constant_time_drift_observation_bias:
        {
            // Test input consistency
            std::shared_ptr< estimatable_parameters::ConstantTimeDriftBiasParameter > biasParameter =
                    std::dynamic_pointer_cast< estimatable_parameters::ConstantTimeDriftBiasParameter >( parameter );
            if( biasParameter == nullptr )
            {
                throw std::runtime_error( "Error, cannot perform bias closure for time drift bias, inconsistent bias types" );
            }

            // Check if bias object is of same type as estimated parameter
            std::shared_ptr< ConstantTimeDriftBias< ObservationSize > > timeBiasObject =
                    std::dynamic_pointer_cast< ConstantTimeDriftBias< ObservationSize > >( observationBias );
            if ( timeBiasObject != nullptr )
            {
                // Check if bias and parameter link properties are equal
                if ( linkEnds == biasParameter->getLinkEnds( ) && observableType == biasParameter->getObservableType( ) )
                {
                    biasParameter->setObservationBiasFunctions(
                            std::bind( &ConstantTimeDriftBias< ObservationSize >::getTemplateFreeConstantObservationBias,
                                       timeBiasObject ),
                            std::bind( &ConstantTimeDriftBias< ObservationSize >::resetConstantObservationBiasTemplateFree,
                                       timeBiasObject, std::placeholders::_1 ) );
                }
            }
            break;
        }
        case estimatable_parameters::arc_wise_time_drift_observation_bias:
        {
            // Test input consistency
            std::shared_ptr< estimatable_parameters::ArcWiseTimeDriftBiasParameter > timeBiasParameter =
                    std::dynamic_pointer_cast< estimatable_parameters::ArcWiseTimeDriftBiasParameter >( parameter );

            if ( timeBiasParameter == nullptr )
            {
                throw std::runtime_error( "Error, cannot perform bias closure for arc-wise time drift biases, inconsistent bias types" );
            }

            // Check if bias object is of same type as estimated parameter
            std::shared_ptr< ArcWiseTimeDriftBias< ObservationSize > > timeBiasObject =
                    std::dynamic_pointer_cast< ArcWiseTimeDriftBias< ObservationSize > >( observationBias );
            if( timeBiasObject != nullptr )
            {
                // Check if bias and parameter link properties are equal
                if ( ( linkEnds == timeBiasParameter->getLinkEnds( ) ) && ( observableType == timeBiasParameter->getObservableType( ) ) &&
                     ( timeBiasObject->getLinkEndIndexForTime( ) == timeBiasParameter->getLinkEndIndex( ) ) &&
                     ( timeBiasObject->getArcStartTimes( ).size( ) == timeBiasParameter->getArcStartTimes( ).size( ) ) )
                {
                    bool doTimesMatch = true;
                    for( unsigned int i = 0 ; i < timeBiasObject->getArcStartTimes( ).size( ) ; i++ )
                    {
                        if ( std::fabs( timeBiasObject->getArcStartTimes( ).at( i ) - timeBiasParameter->getArcStartTimes( ).at( i ) ) >
                             std::max( 1.0E-15 * std::fabs( timeBiasObject->getArcStartTimes( ).at( i ) ),
                                       1.0E-15 * std::fabs( timeBiasParameter->getArcStartTimes( ).at( i ) ) ) )
                        {
                            doTimesMatch = false;
                        }
                    }

                    if( doTimesMatch == true )
                    {
                        timeBiasParameter->setObservationBiasFunctions(
                                std::bind( &ArcWiseTimeDriftBias< ObservationSize >::getTemplateFreeConstantObservationBias,
                                           timeBiasObject ),
                                std::bind( &ArcWiseTimeDriftBias< ObservationSize >::resetConstantObservationBiasTemplateFree,
                                           timeBiasObject, std::placeholders::_1 ) );
                        timeBiasParameter->setLookupScheme( timeBiasObject->getLookupScheme( ) );
                    }
                }
            }
            break;
        }
        case estimatable_parameters::constant_time_observation_bias:
        {
            // Test input consistency
            std::shared_ptr< estimatable_parameters::ConstantTimeBiasParameter > biasParameter =
                    std::dynamic_pointer_cast< estimatable_parameters::ConstantTimeBiasParameter >( parameter );
            if( biasParameter == nullptr )
            {
                throw std::runtime_error( "Error, cannot perform bias closure for time bias, inconsistent bias types" );
            }

            // Check if bias object is of same type as estimated parameter
            std::shared_ptr< ConstantTimeBias< ObservationSize > > timeBiasObject =
                    std::dynamic_pointer_cast< ConstantTimeBias< ObservationSize > >( observationBias );
            if ( timeBiasObject != nullptr )
            {
                // Check if bias and parameter link properties are equal
                if ( linkEnds == biasParameter->getLinkEnds( ) && observableType == biasParameter->getObservableType( ) )
                {
                    biasParameter->setObservationBiasFunctions(
                            std::bind( &ConstantTimeBias< ObservationSize >::getTemplateFreeConstantObservationBias,
                                       timeBiasObject ),
                            std::bind( &ConstantTimeBias< ObservationSize >::resetConstantObservationBiasTemplateFree,
                                       timeBiasObject, std::placeholders::_1 ) );
                }
            }
            break;
        }
        case estimatable_parameters::arc_wise_time_observation_bias:
        {
            // Test input consistency
            std::shared_ptr< estimatable_parameters::ArcWiseTimeBiasParameter > timeBiasParameter =
                    std::dynamic_pointer_cast< estimatable_parameters::ArcWiseTimeBiasParameter >( parameter );

            if ( timeBiasParameter == nullptr )
            {
                throw std::runtime_error( "Error, cannot perform bias closure for arc-wise time biases, inconsistent bias types" );
            }

            // Check if bias object is of same type as estimated parameter
            std::shared_ptr< ArcWiseTimeBias< ObservationSize > > timeBiasObject =
                    std::dynamic_pointer_cast< ArcWiseTimeBias< ObservationSize > >( observationBias );
            if( timeBiasObject != nullptr )
            {
                // Check if bias and parameter link properties are equal
                if ( ( linkEnds == timeBiasParameter->getLinkEnds( ) ) && ( observableType == timeBiasParameter->getObservableType( ) ) &&
                     ( timeBiasObject->getLinkEndIndexForTime( ) == timeBiasParameter->getLinkEndIndex( ) ) &&
                     ( timeBiasObject->getArcStartTimes( ).size( ) == timeBiasParameter->getArcStartTimes( ).size( ) ) )
                {
                    bool doTimesMatch = true;
                    for( unsigned int i = 0 ; i < timeBiasObject->getArcStartTimes( ).size( ) ; i++ )
                    {
                        if ( std::fabs( timeBiasObject->getArcStartTimes( ).at( i ) - timeBiasParameter->getArcStartTimes( ).at( i ) ) >
                             std::max( 1.0E-15 * std::fabs( timeBiasObject->getArcStartTimes( ).at( i ) ),
                                       1.0E-15 * std::fabs( timeBiasParameter->getArcStartTimes( ).at( i ) ) ) )
                        {
                            doTimesMatch = false;
                        }
                    }

                    if( doTimesMatch == true )
                    {
                        timeBiasParameter->setObservationBiasFunctions(
                                std::bind( &ArcWiseTimeBias< ObservationSize >::getTemplateFreeConstantObservationBias,
                                           timeBiasObject ),
                                std::bind( &ArcWiseTimeBias< ObservationSize >::resetConstantObservationBiasTemplateFree,
                                           timeBiasObject, std::placeholders::_1 ) );
                        timeBiasParameter->setLookupScheme( timeBiasObject->getLookupScheme( ) );
                    }
                }
            }
            break;
        }
        default:
            std::string errorMessage = "Error when closing observation bias/estimation loop, did not recognize bias type " +
                    std::to_string( parameter->getParameterName( ).first );
            throw std::runtime_error( errorMessage );

        }
    }
}


template< typename TimeType >
void performTimeBiasPartialClosure(
    const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > timeBiasPartial,
    const std::shared_ptr< propagators::DependentVariablesInterface< TimeType > > dependentVariablesInterface )
{
    std::string bodyName = timeBiasPartial->getParameterName( ).second.first;
    std::shared_ptr< propagators::SingleDependentVariableSaveSettings > totalAccelerationVariable
        = std::make_shared< propagators::SingleDependentVariableSaveSettings >( propagators::total_acceleration_dependent_variable, bodyName );
    std::function< Eigen::VectorXd( const double ) > accelerationPartialFunction =
        [=](const double time)
        {
            return dependentVariablesInterface->getSingleDependentVariable(
                totalAccelerationVariable, time );
        };
    if( std::dynamic_pointer_cast< estimatable_parameters::TimeBiasParameterBase >( timeBiasPartial ) != nullptr )
    {
        std::dynamic_pointer_cast< estimatable_parameters::TimeBiasParameterBase >( timeBiasPartial )->setBodyAccelerationFunction(
            accelerationPartialFunction );
    }
    else
    {
        throw std::runtime_error( "Error when setting time bias parameter closure, type is not supported for parameter " +
            std::to_string( timeBiasPartial->getParameterName( ).first ) + ", " + timeBiasPartial->getParameterName( ).second.first );
    }

}

//! Function to perform the closure between observation models and estimated parameters.
/*!
 *  Function to perform the closure between observation models and estimated parameters. Estimated parameter objects are typically
 *  created prior to observation models. This function must be called for the estimated parameter object creation to be
 *  finalized, in the case link properties are estimated (e.g. observation biases).
 */
template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double >
void performObservationParameterEstimationClosure(
        std::shared_ptr< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType > > observationSimulator ,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > >
        parametersToEstimate,
        std::shared_ptr< propagators::DependentVariablesInterface< TimeType > > dependentVariablesInterface = nullptr )
{
    // Retrieve observation models and parameter
    std::map< LinkEnds, std::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType > > >
            observationModels = observationSimulator->getObservationModels( );
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > > vectorParameters =
            parametersToEstimate->getEstimatedVectorParameters( );

    // Retrieve estimated bias parameters.
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > > vectorBiasParameters;
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > > vectorTimeBiasParameters;

    for( unsigned int i = 0; i < vectorParameters.size( ); i++ )
    {
        if( ( estimatable_parameters::isParameterObservationLinkProperty( vectorParameters.at( i )->getParameterName( ).first ) ) ||
                ( estimatable_parameters::isParameterObservationLinkTimeProperty( vectorParameters.at( i )->getParameterName( ).first ) ) )
        {
            vectorBiasParameters.push_back( vectorParameters.at( i ) );
        }

        if( estimatable_parameters::isParameterObservationLinkTimeProperty( vectorParameters.at( i )->getParameterName( ).first ) )
        {
            vectorTimeBiasParameters.push_back( vectorParameters.at( i ) );
        }
    }

    if( vectorBiasParameters.size( ) > 0 )
    {
        // Retrieve bias objects
        std::map< LinkEnds, std::shared_ptr< ObservationBias< ObservationSize > > > observationBiases =
                extractObservationBiasList( observationModels );

        // Iterate over all combinations of parameters and biases and perform closure for each (if needed)
        for( unsigned int i = 0; i < vectorBiasParameters.size( ); i++ )
        {
            for( typename std::map< LinkEnds, std::shared_ptr< ObservationBias< ObservationSize > > >::const_iterator
                 biasIterator = observationBiases.begin( ); biasIterator != observationBiases.end( ); biasIterator++ )
            {
                performObservationParameterEstimationClosureForSingleModelSet(
                            vectorBiasParameters.at( i ), biasIterator->second, biasIterator->first,
                            observationSimulator->getObservableType( ) );
            }
        }
    }

    for( unsigned int i = 0; i < vectorTimeBiasParameters.size( ); i++ )
    {
        performTimeBiasPartialClosure( vectorTimeBiasParameters.at( i ), dependentVariablesInterface );
    }
}

//! Function to create an object to simulate observations of a given type and associated partials
/*!
 *  Function to create an object to simulate observations of a given type and associated partials
 *  \param observableType Type of observable for which object is to simulate ObservationSimulator
 *  \param observationModelSettingsList Map of settings for the observation models that are to be created in the simulator object: one
 *  for each required set of link ends (each settings object must be consistent with observableType).
 *  \param bodies Map of Body objects that comprise the environment
 *  \param parametersToEstimate Object containing the list of all parameters that are to be estimated
 *  \param stateTransitionMatrixInterface Object used to compute the state transition/sensitivity matrix at a given time
 *  \return Object that simulates the observations of a given type and associated partials
 */
template< int ObservationSize = 1, typename ObservationScalarType, typename TimeType >
std::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType > > createObservationManager(
        const ObservableType observableType,
        const std::vector< std::shared_ptr< ObservationModelSettings  > > observationModelSettingsList,
        const simulation_setup::SystemOfBodies &bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > >
        parametersToEstimate,
        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface >
        stateTransitionMatrixInterface,
        const std::shared_ptr< propagators::DependentVariablesInterface< TimeType > > dependentVariablesInterface =
                std::shared_ptr< propagators::DependentVariablesInterface< TimeType > >( ) )
{
    using namespace observation_models;
    using namespace observation_partials;

    // Create observation simulator
    std::shared_ptr< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType > > observationSimulator =
            createObservationSimulator< ObservationSize, ObservationScalarType, TimeType >(
                observableType, observationModelSettingsList, bodies );

    performObservationParameterEstimationClosure(
                observationSimulator, parametersToEstimate, dependentVariablesInterface );

    // Create observation partials for all link ends/parameters
    std::map< LinkEnds, std::pair< std::map< std::pair< int, int >,
            std::shared_ptr< ObservationPartial< ObservationSize > > >,
            std::shared_ptr< PositionPartialScaling > > > observationPartialsAndScaler;
    if( parametersToEstimate != nullptr )
    {
        observationPartialsAndScaler =
                createObservablePartialsList(
                    observationSimulator->getObservationModels( ), bodies, parametersToEstimate, false, false, dependentVariablesInterface );
    }

    // Split position partial scaling and observation partial objects.
    std::map< LinkEnds, std::map< std::pair< int, int >,
            std::shared_ptr< observation_partials::ObservationPartial< ObservationSize > > > > observationPartials;
    std::map< LinkEnds, std::shared_ptr< observation_partials::PositionPartialScaling  > > observationPartialScalers;
    splitObservationPartialsAndScalers( observationPartialsAndScaler, observationPartials, observationPartialScalers );


    return std::make_shared< ObservationManager< ObservationSize, ObservationScalarType, TimeType > >(
                observableType, observationSimulator, observationPartials,
                observationPartialScalers, stateTransitionMatrixInterface, dependentVariablesInterface );
}


//! Function to create an object to simulate observations of a given type and associated partials
/*!
 *  Function to create an object to simulate observations of a given type and associated partials
 *  \param observableType Type of observable for which object is to simulate ObservationSimulator
 *  \param observationModelSettingsList Map of settings for the observation models that are to be created in the simulator object: one
 *  for each required set of link ends (each settings object must be consistent with observableType).
 *  \param bodies Map of Body objects that comprise the environment
 *  \param parametersToEstimate Object containing the list of all parameters that are to be estimated
 *  \param stateTransitionMatrixInterface Object used to compute the state transition/sensitivity matrix at a given time
 *  \return Object that simulates the observations of a given type and associated partials
 */
template< typename ObservationScalarType, typename TimeType >
std::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType > > createObservationManagerBase(
        const ObservableType observableType,
        const std::vector< std::shared_ptr< ObservationModelSettings > > observationModelSettingsList,
        const simulation_setup::SystemOfBodies &bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > parametersToEstimate,
        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface,
        const std::shared_ptr< propagators::DependentVariablesInterface< TimeType > > dependentVariablesInterface =
                std::shared_ptr< propagators::DependentVariablesInterface< TimeType > >( ) )
{
    std::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType > > observationManager;
    switch( observableType )
    {
    case one_way_range:
    case n_way_range:
    case one_way_doppler:
    case two_way_doppler:
    case one_way_differenced_range:
    case n_way_differenced_range:
    case dsn_one_way_averaged_doppler:
    case dsn_n_way_averaged_doppler:
        observationManager = createObservationManager< 1, ObservationScalarType, TimeType >(
                    observableType, observationModelSettingsList, bodies, parametersToEstimate,
                    stateTransitionMatrixInterface, dependentVariablesInterface );
        break;
    case angular_position:
    case relative_angular_position:
        observationManager = createObservationManager< 2, ObservationScalarType, TimeType >(
                    observableType, observationModelSettingsList, bodies, parametersToEstimate,
                    stateTransitionMatrixInterface, dependentVariablesInterface );
        break;
    case position_observable:
    case euler_angle_313_observable:
    case velocity_observable:
        observationManager = createObservationManager< 3, ObservationScalarType, TimeType >(
                    observableType, observationModelSettingsList, bodies, parametersToEstimate,
                    stateTransitionMatrixInterface, dependentVariablesInterface );
        break;
    case relative_position_observable:
        observationManager = createObservationManager< 3, ObservationScalarType, TimeType >(
                observableType, observationModelSettingsList, bodies, parametersToEstimate,
                stateTransitionMatrixInterface, dependentVariablesInterface );
        break;
    default:
        throw std::runtime_error(
                    "Error when making observation manager, could not identify observable type " +
                    std::to_string( observableType ) );
    }
    return observationManager;
}

template< typename ObservationScalarType, typename TimeType >
std::map< ObservableType, std::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType > > > createObservationManagersBase(
    const std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >& observationSettingsList,
    const simulation_setup::SystemOfBodies &bodies,
    const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > fullParameters,
    const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface,
    const std::shared_ptr< propagators::DependentVariablesInterface< TimeType > > dependentVariablesInterface =
    std::shared_ptr< propagators::DependentVariablesInterface< TimeType > >( ) )
{
    std::map< ObservableType, std::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType > > > observationManagers;
    std::map< ObservableType, std::vector< std::shared_ptr< ObservationModelSettings > > > sortedObservationSettingsList =
        sortObservationModelSettingsByType( observationSettingsList );

    for( auto it : sortedObservationSettingsList )
    {
        // Call createObservationSimulator of required observation size
        ObservableType observableType = it.first;

        // Create observation manager for current observable.
        observationManagers[ observableType ] = createObservationManagerBase< ObservationScalarType, TimeType >(
            observableType, it.second, bodies, fullParameters,
            stateTransitionMatrixInterface, dependentVariablesInterface );
    }
    return observationManagers;
}

}


}


#endif // TUDAT_CREATEOBSERVATIONMANAGER_H

/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEOBSERVATIONBIASPARTIALS_H
#define TUDAT_CREATEOBSERVATIONBIASPARTIALS_H

#include <vector>
#include <map>

#include <memory>

#include <Eigen/Core>

#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/observationBiasPartial.h"

namespace tudat
{

namespace observation_partials
{


template< int ObservationSize >
std::shared_ptr< ObservationPartial< ObservationSize > > getPartialWrtBodyTranslationalStateFromList(
    const std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > >& observationPartials,
    const std::string& bodyName )
{
    std::shared_ptr< ObservationPartial< ObservationSize > > selectedObservationPartial = nullptr;
    for( auto it : observationPartials )
    {
        estimatable_parameters::EstimatebleParameterIdentifier currentParameter = it.second->getParameterIdentifier( );
        if( currentParameter.first == estimatable_parameters::initial_body_state && currentParameter.second.second == bodyName )
        {
            if( selectedObservationPartial != nullptr )
            {
                throw std::runtime_error( "Error when getting observation partial w.r.t. translational state of " + bodyName + ", found several options" );
            }
            selectedObservationPartial = it.second;
        }
    }
    return selectedObservationPartial;
}

template< int ObservationSize >
std::shared_ptr< ObservationPartial< ObservationSize > > getPartialWrtBodyTranslationalState(
    const std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > >& observationPartials,
    const std::string& bodyName,
    std::function< std::shared_ptr< ObservationPartial< ObservationSize > >( const std::string& ) > partialCreationFunction )
{
    std::shared_ptr< ObservationPartial< ObservationSize > > observationPartial = getPartialWrtBodyTranslationalStateFromList(
        observationPartials, bodyName );
    if( observationPartial == nullptr && partialCreationFunction != nullptr )
    {
        observationPartial = partialCreationFunction( bodyName );
    }
    else if( observationPartial == nullptr && partialCreationFunction == nullptr )
    {
        throw std::runtime_error( "Error when getting partial of observable w.r.t. state of body " + bodyName + " for bias partial." );
    }
    return observationPartial;
}
//! Function to create partials of observation w.r.t. a link property.
/*!
 *  Function to create partials of observation w.r.t. a link property, e.g. a parameter that does not influence either link end's
 *  dynamics, only the observable itself, such as observation biases and clock parameters.
 *  \param linkEnds Link ends of observable for which partial is to be made.
 *  \param observableType Type of observable for which partial is to be made.
 *  \param parameterToEstimate Parameter w.r.t. which the partial is to be taken
 *  \param useBiasPartials Boolean to denote whether this function should create partials w.r.t. observation bias parameters
 *  \return Object that computes the partial of the observation w.r.t. parameterToEstimate (nullptr if no dependency).
 */
template< int ObservationSize >
std::shared_ptr< ObservationPartial< ObservationSize > > createObservationPartialWrtLinkProperty(
    const observation_models::LinkEnds& linkEnds,
    const observation_models::ObservableType observableType,
    const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameterToEstimate,
    const bool useBiasPartials = true,
    const std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > >& observationPartials =
        std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > >( ),
    const std::function< std::shared_ptr< ObservationPartial< ObservationSize > >( const std::string& ) > partialWrtStateCreationFunction = nullptr )
{
//    std::cout << "observableType: " << observableType << "\n\n";

    std::shared_ptr< ObservationPartial< ObservationSize > > observationPartial;

    // Check parameter type
    switch( parameterToEstimate->getParameterName( ).first )
    {
    case estimatable_parameters::constant_additive_observation_bias:
    {
        if( useBiasPartials )
        {
            // Check input consistency
            std::shared_ptr< estimatable_parameters::ConstantObservationBiasParameter > constantBias =
                std::dynamic_pointer_cast< estimatable_parameters::ConstantObservationBiasParameter >(
                    parameterToEstimate );
            if( constantBias == nullptr )
            {
                throw std::runtime_error( "Error when making partial w.r.t. observation bias, type is inconsistent" );
            }
            else
            {
                // Check dependency between parameter and link properties.
                if( linkEnds == constantBias->getLinkEnds( ) && observableType == constantBias->getObservableType( ) )
                {
                    observationPartial = std::make_shared< ObservationPartialWrtConstantAbsoluteBias< ObservationSize > >(
                        observableType, linkEnds );
                }
            }
        }
        break;
    }
    case estimatable_parameters::arcwise_constant_additive_observation_bias:
    {
        if( useBiasPartials )
        {
            // Check input consistency
            std::shared_ptr< estimatable_parameters::ArcWiseObservationBiasParameter > arcwiseBias =
                std::dynamic_pointer_cast< estimatable_parameters::ArcWiseObservationBiasParameter >(
                    parameterToEstimate );
            if( arcwiseBias == nullptr )
            {
                throw std::runtime_error( "Error when making partial w.r.t. arcwise observation bias, type is inconsistent" );
            }
            else
            {
                // Check dependency between parameter and link properties.
                if( linkEnds == arcwiseBias->getLinkEnds( ) && observableType == arcwiseBias->getObservableType( ) )
                {
                    observationPartial = std::make_shared< ObservationPartialWrtArcWiseAbsoluteBias< ObservationSize > >(
                        observableType, linkEnds,
                        arcwiseBias->getLookupScheme( ),
                        arcwiseBias->getLinkEndIndex( ),
                        arcwiseBias->getArcStartTimes( ).size( ) );
                }
            }
        }
        break;
    }
    case estimatable_parameters::constant_relative_observation_bias:
    {
        if( useBiasPartials )
        {
            // Check input consistency
            std::shared_ptr< estimatable_parameters::ConstantObservationBiasParameter > constantBias =
                std::dynamic_pointer_cast< estimatable_parameters::ConstantObservationBiasParameter >(
                    parameterToEstimate );
            if( constantBias == nullptr )
            {
                throw std::runtime_error( "Error when making partial w.r.t. observation bias, type is inconsistent" );
            }
            else
            {
                // Check dependency between parameter and link properties.
                if( linkEnds == constantBias->getLinkEnds( ) && observableType == constantBias->getObservableType( ) )
                {
                    observationPartial = std::make_shared< ObservationPartialWrtConstantRelativeBias< ObservationSize > >(
                        observableType, linkEnds );

                }
            }
        }
        break;
    }
    case estimatable_parameters::arcwise_constant_relative_observation_bias:
    {
        if( useBiasPartials )
        {
            // Check input consistency
            std::shared_ptr< estimatable_parameters::ArcWiseObservationBiasParameter > arcwiseBias =
                std::dynamic_pointer_cast< estimatable_parameters::ArcWiseObservationBiasParameter >(
                    parameterToEstimate );
            if( arcwiseBias == nullptr )
            {
                throw std::runtime_error( "Error when making partial w.r.t. arcwise relative observation bias, type is inconsistent" );
            }
            else
            {
                // Check dependency between parameter and link properties.
                if( linkEnds == arcwiseBias->getLinkEnds( ) && observableType == arcwiseBias->getObservableType( ) )
                {
                    observationPartial = std::make_shared< ObservationPartialWrtArcWiseRelativeBias< ObservationSize > >(
                        observableType, linkEnds,
                        arcwiseBias->getLookupScheme( ),
                        arcwiseBias->getLinkEndIndex( ),
                        arcwiseBias->getArcStartTimes( ).size( ) );
                }
            }
        }
        break;
    }
//    case estimatable_parameters::constant_time_drift_observation_bias:
//    {
//        if( useBiasPartials )
//        {
//            // Check input consistency
//            std::shared_ptr< estimatable_parameters::ConstantTimeDriftBiasParameter > constantTimeDriftBias =
//                std::dynamic_pointer_cast< estimatable_parameters::ConstantTimeDriftBiasParameter >(
//                    parameterToEstimate );
//            if( constantTimeDriftBias == nullptr )
//            {
//                throw std::runtime_error( "Error when making partial w.r.t. time drift bias, type is inconsistent" );
//            }
//            else
//            {
//                // Check dependency between parameter and link properties.
//                if( linkEnds == constantTimeDriftBias->getLinkEnds( ) && observableType == constantTimeDriftBias->getObservableType( ) )
//                {
//                    observationPartial = std::make_shared< ObservationPartialWrtConstantTimeDriftBias< ObservationSize > >(
//                        observableType, linkEnds, constantTimeDriftBias->getLinkEndIndex( ), constantTimeDriftBias->getReferenceEpoch( ) );
//                }
//            }
//        }
//        break;
//    }
//    case estimatable_parameters::arc_wise_time_drift_observation_bias:
//    {
//        if( useBiasPartials )
//        {
//            // Check input consistency
//            std::shared_ptr< estimatable_parameters::ArcWiseTimeDriftBiasParameter > arcwiseTimeDriftBias =
//                std::dynamic_pointer_cast< estimatable_parameters::ArcWiseTimeDriftBiasParameter >( parameterToEstimate );
//            if ( arcwiseTimeDriftBias == nullptr )
//            {
//                throw std::runtime_error( "Error when making partial w.r.t. arcwise time drift bias, type is inconsistent" );
//            }
//            else
//            {
//                // Check dependency between parameter and link properties.
//                if ( linkEnds == arcwiseTimeDriftBias->getLinkEnds( ) && observableType == arcwiseTimeDriftBias->getObservableType( ) )
//                {
//                    observationPartial = std::make_shared< ObservationPartialWrtArcWiseTimeDriftBias< ObservationSize > >(
//                        observableType, linkEnds,
//                        arcwiseTimeDriftBias->getLookupScheme( ),
//                        arcwiseTimeDriftBias->getLinkEndIndex( ),
//                        arcwiseTimeDriftBias->getArcStartTimes( ).size( ),
//                        arcwiseTimeDriftBias->getReferenceEpochs( ) );
//                }
//            }
//        }
//        break;
//    }
    case estimatable_parameters::constant_time_observation_bias:
    {
        if( useBiasPartials )
        {
            // Check input consistency
            std::shared_ptr< estimatable_parameters::ConstantTimeBiasParameter > constantTimeBias =
                std::dynamic_pointer_cast< estimatable_parameters::ConstantTimeBiasParameter >( parameterToEstimate );
            if( constantTimeBias == nullptr )
            {
                throw std::runtime_error( "Error when making partial w.r.t. time bias, type is inconsistent" );
            }
            else
            {
                // Check dependency between parameter and link properties.
                if( linkEnds == constantTimeBias->getLinkEnds( ) && observableType == constantTimeBias->getObservableType( ) )
                {
                    std::shared_ptr< ObservationPartial< ObservationSize > > observationWrtTransmitterStatePartial = getPartialWrtBodyTranslationalState(
                        observationPartials, linkEnds.at( observation_models::transmitter ).bodyName_, partialWrtStateCreationFunction );
                    std::shared_ptr< ObservationPartial< ObservationSize > > observationWrtReceiverStatePartial = getPartialWrtBodyTranslationalState(
                        observationPartials, linkEnds.at( observation_models::receiver ).bodyName_, partialWrtStateCreationFunction );


                    std::shared_ptr< TimeBiasPartial< ObservationSize > > timeBiasPartial= std::make_shared< TimeBiasPartial< ObservationSize > >(
                       linkEnds, observableType,
                       constantTimeBias->getLinkEndId( ),
                       constantTimeBias->getReferenceLinkEnd( ),
                       observationWrtTransmitterStatePartial,
                       observationWrtReceiverStatePartial,
                       constantTimeBias->getBodyAccelerationFunction( ) );
//                    partialWrtParameterBodyState
                    observationPartial = std::make_shared< ObservationPartialWrtConstantTimeBias< ObservationSize > >( timeBiasPartial );
                }
            }
        }
        break;
    }
//    case estimatable_parameters::arc_wise_time_observation_bias:
//    {
//        if( useBiasPartials )
//        {
//            // Check input consistency
//            std::shared_ptr< estimatable_parameters::ArcWiseTimeBiasParameter > arcwiseTimeBias =
//                std::dynamic_pointer_cast< estimatable_parameters::ArcWiseTimeBiasParameter >( parameterToEstimate );
//            if ( arcwiseTimeBias == nullptr )
//            {
//                throw std::runtime_error( "Error when making partial w.r.t. arcwise time bias, type is inconsistent" );
//            }
//            else
//            {
//                // Check dependency between parameter and link properties.
//                if ( linkEnds == arcwiseTimeBias->getLinkEnds( ) && observableType == arcwiseTimeBias->getObservableType( ) )
//                {
//                    observationPartial = std::make_shared< ObservationPartialWrtArcWiseTimeBias< ObservationSize > >(
//                        observableType, linkEnds,
//                        arcwiseTimeBias->getLookupScheme( ),
//                        arcwiseTimeBias->getLinkEndIndex( ),
//                        arcwiseTimeBias->getArcStartTimes( ).size( ) );
//                }
//            }
//        }
//        break;
//    }
    default:
        throw std::runtime_error( "Error when making partial w.r.t. link property, parameter " + std::to_string( parameterToEstimate->getParameterName( ).first ) + " not found." );
    }

    return observationPartial;
}

}

}


#endif // TUDAT_CREATEOBSERVATIONBIASPARTIALS_H

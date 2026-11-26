/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SETINITIALSTATEVECTORFROMPARAMETERSET_H
#define TUDAT_SETINITIALSTATEVECTORFROMPARAMETERSET_H

#include <algorithm>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameterSet.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to get initial state vector of estimated dynamical states.
/*!
 *  Function to get initial state vector of estimated dynamical states (i.e. presently estimated state at propagation
 *  start time.
 *  \param estimatableParameters Object containing all parameters that are to be estimated.
 *  \param propagatorSettings Object containing propagation settings to be used
 *  \return State vector of estimated dynamics at propagation start time.
 */
template< typename InitialStateParameterType = double, typename TimeType = double >
void setInitialStateVectorFromParameterSet(
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > > estimatableParameters,
        const std::shared_ptr< propagators::PropagatorSettings< InitialStateParameterType > > propagatorSettings )
{
    typedef Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > VectorType;

    // Retrieve initial dynamical parameters.
    std::vector< std::shared_ptr<
            estimatable_parameters::EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
            initialDynamicalParameters = estimatableParameters->getEstimatedInitialStateParameters( );

    // Initialize state vector.
    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialStateVector =
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >::Zero(
                    estimatableParameters->getInitialDynamicalStateParameterSize( ), 1 );

    if( std::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< InitialStateParameterType, TimeType > >(
                propagatorSettings ) != nullptr )
    {
        int vectorSize = 0;
        // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
        for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
        {
            if( isParameterDynamicalPropertyInitialState( initialDynamicalParameters.at( i )->getParameterName( ).first ) )
            {
                int currentParameterSize = initialDynamicalParameters.at( i )->getParameterSize( );
                initialStateVector.block( vectorSize, 0, currentParameterSize, 1 ) =
                        initialDynamicalParameters.at( i )->getParameterValue( );

                vectorSize += currentParameterSize;
            }
        }

        propagatorSettings->resetInitialStates( initialStateVector.block( 0, 0, vectorSize, 1 ) );
    }
    else if( std::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< InitialStateParameterType, TimeType > >(
                     propagatorSettings ) )
    {
        std::shared_ptr< propagators::MultiArcPropagatorSettings< InitialStateParameterType, TimeType > > multiArcSettings =
                std::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< InitialStateParameterType, TimeType > >(
                        propagatorSettings );

        std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< InitialStateParameterType, TimeType > > >
                singleArcSettings = multiArcSettings->getSingleArcSettings( );
        int numberOfArcs = singleArcSettings.size( );

        // Counting in how many arcs each body has already been propagated/estimated.
        std::map< propagators::IntegratedStateType, std::map< std::string, unsigned int > > initialStatesBodiesCounter;

        for( int i = 0; i < numberOfArcs; i++ )
        {
            std::map< propagators::IntegratedStateType, std::map< std::pair< std::string, std::string >, VectorType > >
                    currentArcInitialStates;

            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< InitialStateParameterType, TimeType > >
                    singleArcTranslationalStatePropagatorSettings = std::dynamic_pointer_cast<
                            propagators::TranslationalStatePropagatorSettings< InitialStateParameterType, TimeType > >(
                            singleArcSettings.at( i ) );

            if( singleArcTranslationalStatePropagatorSettings == nullptr )
            {
                throw std::runtime_error( "Propagator settings for arc " + std::to_string( i ) +
                                          " are not translational state propagator settings. "
                                          " Other propagator settings not supported (yet) for multi-arc propagation." );
            }
            else
            {
                for( unsigned int j = 0; j < initialDynamicalParameters.size( ); j++ )
                {
                    VectorType currentParameterValue = initialDynamicalParameters.at( j )->getParameterValue( );
                    int currentParameterSize = initialDynamicalParameters.at( j )->getParameterSize( );
                    std::pair< std::string, std::string > bodyIdentifier = initialDynamicalParameters.at( j )->getParameterName( ).second;
                    //                    std::cout << "body identifier: " << bodyIdentifier.first << " & " << bodyIdentifier.second <<
                    //                    "\n\n";

                    auto itr = std::find( singleArcTranslationalStatePropagatorSettings->bodiesToIntegrate_.begin( ),
                                          singleArcTranslationalStatePropagatorSettings->bodiesToIntegrate_.end( ),
                                          bodyIdentifier.first );
                    if( itr != singleArcTranslationalStatePropagatorSettings->bodiesToIntegrate_.cend( ) )
                    {
                        //                        std::cout << "body " << bodyIdentifier.first << " - detected for arc " << i << "\n\n";

                        switch( initialDynamicalParameters.at( j )->getParameterName( ).first )
                        {
                            case estimatable_parameters::arc_wise_initial_body_state: {
                                //                                std::cout << "current parameter size: " << currentParameterSize << "\n\n";

                                if( ( initialStatesBodiesCounter.count( propagators::translational_state ) == 0 ) ||
                                    initialStatesBodiesCounter.at( propagators::translational_state ).count( bodyIdentifier.first ) == 0 )
                                {
                                    initialStatesBodiesCounter[ propagators::translational_state ][ bodyIdentifier.first ] = 0;
                                }
                                int index = initialStatesBodiesCounter.at( propagators::translational_state ).at( bodyIdentifier.first );
                                //                                std::cout << "index: " << index << "\n\n";

                                //                                if ( currentParameterSize / numberOfArcs != 6 )
                                //                                {
                                //                                    throw std::runtime_error("Error when moving initial states from
                                //                                    parameters to propagator settings. Incompatible multi-arc
                                //                                    translational state size found");
                                //                                }
                                //                                std::cout << "full parameters values: " <<
                                //                                currentParameterValue.transpose( ) << "\n\n"; std::cout << "subset
                                //                                parameters values: " << currentParameterValue.segment( index * 6, 6
                                //                                ).transpose( ) << "\n\n";
                                currentArcInitialStates[ propagators::translational_state ][ bodyIdentifier ] =
                                        currentParameterValue.segment( index * 6, 6 );

                                // update counter of propagated/estimated bodies
                                initialStatesBodiesCounter.at( propagators::translational_state ).at( bodyIdentifier.first )++;

                                break;
                            }
                            default:
                                throw std::runtime_error(
                                        "Error when moving initial states from parameters to propagator settings. Multi-arc parameter type "
                                        "not recognized" );
                        }
                    }
                }
            }
            propagators::resetSingleArcInitialStates( singleArcSettings.at( i ), currentArcInitialStates );
            multiArcSettings->updateInitialState( );
        }
    }
    else if( std::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< InitialStateParameterType, TimeType > >(
                     propagatorSettings ) )
    {
        std::shared_ptr< propagators::HybridArcPropagatorSettings< InitialStateParameterType, TimeType > > hybridArcSettings =
                std::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< InitialStateParameterType, TimeType > >(
                        propagatorSettings );

        std::shared_ptr< propagators::SingleArcPropagatorSettings< InitialStateParameterType, TimeType > > singleArcSettings =
                hybridArcSettings->getSingleArcPropagatorSettings( );
        std::shared_ptr< propagators::MultiArcPropagatorSettings< InitialStateParameterType, TimeType > > multiArcSettings =
                hybridArcSettings->getMultiArcPropagatorSettings( );

        setInitialStateVectorFromParameterSet< InitialStateParameterType, TimeType >(
                std::make_shared< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >(
                        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > >( ),
                        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > >( ),
                        estimatableParameters->getEstimatedSingleArcInitialStateParameters( ) ),
                singleArcSettings );
        setInitialStateVectorFromParameterSet< InitialStateParameterType, TimeType >(
                std::make_shared< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >(
                        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > >( ),
                        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > >( ),
                        estimatableParameters->getEstimatedMultiArcInitialStateParameters( ) ),
                multiArcSettings );

        hybridArcSettings->updateInitialState( );
    }
    else
    {
        throw std::runtime_error( "Error when setting initial state vector from parameter set; type not identified" );
    }
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_SETINITIALSTATEVECTORFROMPARAMETERSET_H

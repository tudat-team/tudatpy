/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#ifndef CREATEOBSERVATIONMANAGER_H
#define CREATEOBSERVATIONMANAGER_H

#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/ObservationModels/observationManager.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/observationBiasParameter.h"
#include "Tudat/SimulationSetup/EstimationSetup/createObservationModel.h"
#include "Tudat/SimulationSetup/EstimationSetup/createObservationPartials.h"
#include "Tudat/Astrodynamics/ObservationModels/oneWayRangeObservationModel.h"

namespace tudat
{

namespace observation_models
{

//! Function to create an object to simulate observations of a given type
/*!
 *  Function to create an object to simulate observations of a given type
 *  \param observableType Type of observable for which object is to simulate ObservationSimulator
 *  \param settingsPerLinkEnds Map of settings for the observation models that are to be created in the simulator object: one
 *  for each required set of link ends (each settings object must be consistent with observableType).
 *  \param bodyMap Map of Body objects that comprise the environment
 *  \return Object that simulates the observables according to the provided settings.
 */
template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double >
boost::shared_ptr< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType > > createObservationSimulator(
        const ObservableType observableType,
        const std::map< LinkEnds, boost::shared_ptr< ObservationSettings  > > settingsPerLinkEnds,
        const simulation_setup::NamedBodyMap &bodyMap )
{
    std::map< LinkEnds, boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType > > >
            observationModels;

    // Iterate over all link ends
    for( std::map< LinkEnds, boost::shared_ptr< ObservationSettings  > >::const_iterator settingIterator =
         settingsPerLinkEnds.begin( ); settingIterator != settingsPerLinkEnds.end( ); settingIterator++ )
    {
        observationModels[ settingIterator->first ] = ObservationModelCreator<
                ObservationSize, ObservationScalarType, TimeType >::createObservationModel(
                    settingIterator->first, settingIterator->second, bodyMap );
    }

    return boost::make_shared< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType > >(
                observableType, observationModels );


}



template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double >
void performObservationParameterEstimationClosure(
        boost::shared_ptr< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType > > observationSimulator ,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > >
        parametersToEstimate )
{
    std::map< LinkEnds, boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType > > > observationModels =
            observationSimulator->getObservationModels( );

    for( unsigned int i = 0; i < parametersToEstimate->getEstimatedVectorParameters( ).size( ); i++ )
    {
        if( parametersToEstimate->getEstimatedVectorParameters( ).at( i )->getParameterName( ).first ==
                estimatable_parameters::constant_additive_observation_bias )
        {
            boost::shared_ptr< estimatable_parameters::ConstantObservationBiasParameter > biasParameter =
                    boost::dynamic_pointer_cast< estimatable_parameters::ConstantObservationBiasParameter >(
                        parametersToEstimate->getEstimatedVectorParameters( ).at( i ) );

            if( observationSimulator->getObservableType( ) == biasParameter->getObservableType( ) )
            {
                if( observationModels.count( biasParameter->getLinkEnds( ) ) != 0 )
                {
                    if( observationModels.at( biasParameter->getLinkEnds( ) )->getObservationBiasCalculator( ) != NULL )
                    {
                        boost::shared_ptr< ConstantObservationBias< ObservationSize > > constantBiasObject =
                                boost::dynamic_pointer_cast< ConstantObservationBias< ObservationSize > >(
                                    observationModels.at( biasParameter->getLinkEnds( ) )->getObservationBiasCalculator( ) );
                        if( constantBiasObject != NULL )
                        {
                            biasParameter->setObservationBiasFunctions(
                                        boost::bind( &ConstantObservationBias< ObservationSize >::getTemplateFreeConstantObservationBias,
                                                     constantBiasObject ),
                                        boost::bind( &ConstantObservationBias< ObservationSize >::resetConstantObservationBiasTemplateFree,
                                                     constantBiasObject, _1 ) );
                        }
                        else
                        {
                            throw std::runtime_error( "Error, cannot perform bias closure for additive bias, inconsistent bias types" );
                        }
                    }
                    else
                    {
                        throw std::runtime_error( "Error, cannot perform bias closure for additive bias" );
                    }
                }
            }
        }
        else if( parametersToEstimate->getEstimatedVectorParameters( ).at( i )->getParameterName( ).first ==
                 estimatable_parameters::constant_relative_observation_bias )
        {
            boost::shared_ptr< estimatable_parameters::ConstantRelativeObservationBiasParameter > biasParameter =
                    boost::dynamic_pointer_cast< estimatable_parameters::ConstantRelativeObservationBiasParameter >(
                        parametersToEstimate->getEstimatedVectorParameters( ).at( i ) );

            if( observationSimulator->getObservableType( ) == biasParameter->getObservableType( ) )
            {
                if( observationModels.count( biasParameter->getLinkEnds( ) ) != 0 )
                {
                    if( observationModels.at( biasParameter->getLinkEnds( ) )->getObservationBiasCalculator( ) != NULL )
                    {
                        boost::shared_ptr< ConstantRelativeObservationBias< ObservationSize > > constantBiasObject =
                                boost::dynamic_pointer_cast< ConstantRelativeObservationBias< ObservationSize > >(
                                    observationModels.at( biasParameter->getLinkEnds( ) )->getObservationBiasCalculator( ) );
                        if( constantBiasObject != NULL )
                        {
                            biasParameter->setObservationBiasFunctions(
                                        boost::bind( &ConstantRelativeObservationBias< ObservationSize >::getTemplateFreeConstantObservationBias,
                                                     constantBiasObject ),
                                        boost::bind( &ConstantRelativeObservationBias< ObservationSize >::resetConstantObservationBiasTemplateFree,
                                                     constantBiasObject, _1 ) );                        }
                        else
                        {
                            throw std::runtime_error( "Error, cannot perform bias closure for additive bias, inconsistent bias types" );
                        }
                    }
                    else
                    {
                        throw std::runtime_error( "Error, cannot perform bias closure for additive bias" );
                    }
                }
            }

        }
    }

}

//! Function to create an object to simulate observations of a given type and associated partials
/*!
 *  Function to create an object to simulate observations of a given type and associated partials
 *  \param observableType Type of observable for which object is to simulate ObservationSimulator
 *  \param settingsPerLinkEnds Map of settings for the observation models that are to be created in the simulator object: one
 *  for each required set of link ends (each settings object must be consistent with observableType).
 *  \param bodyMap Map of Body objects that comprise the environment
 *  \param parametersToEstimate Object containing the list of all parameters that are to be estimated
 *  \param stateTransitionMatrixInterface Object used to compute the state transition/sensitivity matrix at a given time
 *  \return Object that simulates the observations of a given type and associated partials
 */
template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double >
boost::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType > > createObservationManager(
        const ObservableType observableType,
        const std::map< LinkEnds, boost::shared_ptr< ObservationSettings  > > settingsPerLinkEnds,
        const simulation_setup::NamedBodyMap &bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > >
        parametersToEstimate,
        const boost::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface >
        stateTransitionMatrixInterface )
{
    using namespace observation_models;
    using namespace observation_partials;

    // Create observation simulator
    boost::shared_ptr< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType > > observationSimulator =
            createObservationSimulator< ObservationSize, ObservationScalarType, TimeType >(
                observableType, settingsPerLinkEnds, bodyMap );

    performObservationParameterEstimationClosure(
                observationSimulator, parametersToEstimate );

    // Get light-time corrections for current observable
    PerLinkEndPerLightTimeSolutionCorrections lightTimeCorrectionList =
            getLightTimeCorrectionsList( observationSimulator->getObservationModels( ) );

    // Create observation partials for all link ends/parameters
    boost::shared_ptr< ObservationPartialCreator< ObservationSize, ObservationScalarType > > observationPartialCreator;
    std::map< LinkEnds, std::pair< std::map< std::pair< int, int >,
            boost::shared_ptr< ObservationPartial< ObservationSize > > >,
            boost::shared_ptr< PositionPartialScaling > > > observationPartialsAndScaler;
    if( parametersToEstimate != NULL )
    {
        observationPartialsAndScaler =
                observationPartialCreator->createObservationPartials(
                    observableType, utilities::createVectorFromMapKeys( settingsPerLinkEnds ), bodyMap, parametersToEstimate,
                    lightTimeCorrectionList );
    }

    // Split position partial scaling and observation partial objects.
    std::map< LinkEnds, std::map< std::pair< int, int >,
            boost::shared_ptr< observation_partials::ObservationPartial< ObservationSize > > > > observationPartials;
    std::map< LinkEnds, boost::shared_ptr< observation_partials::PositionPartialScaling  > > observationPartialScalers;
    splitObservationPartialsAndScalers( observationPartialsAndScaler, observationPartials, observationPartialScalers );


    return boost::make_shared< ObservationManager< ObservationSize, ObservationScalarType, TimeType > >(
                observableType, observationSimulator, observationPartials,
                observationPartialScalers, stateTransitionMatrixInterface );
}


//! Function to create an object to simulate observations of a given type and associated partials
/*!
 *  Function to create an object to simulate observations of a given type and associated partials
 *  \param observableType Type of observable for which object is to simulate ObservationSimulator
 *  \param settingsPerLinkEnds Map of settings for the observation models that are to be created in the simulator object: one
 *  for each required set of link ends (each settings object must be consistent with observableType).
 *  \param bodyMap Map of Body objects that comprise the environment
 *  \param parametersToEstimate Object containing the list of all parameters that are to be estimated
 *  \param stateTransitionMatrixInterface Object used to compute the state transition/sensitivity matrix at a given time
 *  \return Object that simulates the observations of a given type and associated partials
 */
template< typename ObservationScalarType = double, typename TimeType = double >
boost::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType > > createObservationManagerBase(
        const ObservableType observableType,
        const std::map< LinkEnds, boost::shared_ptr< ObservationSettings  > > settingsPerLinkEnds,
        const simulation_setup::NamedBodyMap &bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > parametersToEstimate,
        const boost::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface )
{
    boost::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType > > observationManager;
    switch( observableType )
    {
    case one_way_range:
        observationManager = createObservationManager< 1, ObservationScalarType, TimeType >(
                    observableType, settingsPerLinkEnds, bodyMap, parametersToEstimate,
                    stateTransitionMatrixInterface );
        break;
    case one_way_doppler:
        observationManager = createObservationManager< 1, ObservationScalarType, TimeType >(
                    observableType, settingsPerLinkEnds, bodyMap, parametersToEstimate,
                    stateTransitionMatrixInterface );
        break;
    case one_way_differenced_range:
        observationManager = createObservationManager< 1, ObservationScalarType, TimeType >(
                    observableType, settingsPerLinkEnds, bodyMap, parametersToEstimate,
                    stateTransitionMatrixInterface );
        break;
    case angular_position:
        observationManager = createObservationManager< 2, ObservationScalarType, TimeType >(
                    observableType, settingsPerLinkEnds, bodyMap, parametersToEstimate,
                    stateTransitionMatrixInterface );
        break;
    case position_observable:
        observationManager = createObservationManager< 3, ObservationScalarType, TimeType >(
                    observableType, settingsPerLinkEnds, bodyMap, parametersToEstimate,
                    stateTransitionMatrixInterface );
        break;
    default:
        throw std::runtime_error(
                    "Error when making observation manager, could not identify observable type " +
                    boost::lexical_cast< std::string >( observableType ) );
    }
    return observationManager;
}

}


}


#endif // CREATEOBSERVATIONMANAGER_H

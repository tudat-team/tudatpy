/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATION_COLLECTION_CONVERSION_H
#define TUDAT_OBSERVATION_COLLECTION_CONVERSION_H

#include <memory>
#include <stdexcept>
#include <vector>

#include "tudat/simulation/estimation_setup/observationCollection.h"
#include "tudat/simulation/estimation_setup/observationDataset.h"

/*!
 * Conversion between the deprecated SingleObservationSet/ObservationCollection classes and ObservationDataset.
 * All observation functionality is implemented for ObservationDataset. Functions that retain an ObservationCollection
 * interface convert the collection to a dataset, call the dataset implementation, and write the results back into the
 * collection (or convert a resulting dataset to a new collection).
 */

namespace tudat
{

namespace observation_models
{

//! Add (a copy of) a deprecated single observation set to an observation dataset, returning the id of the new dataset set.
template< typename ObservationScalarType, typename TimeType >
int addSingleObservationSetToDataset( ObservationDataset< ObservationScalarType, TimeType >& observationDataset,
                                      const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > >& observationSet )
{
    if( observationSet == nullptr )
    {
        throw std::runtime_error( "Error when adding single observation set to observation dataset, input set is null." );
    }
    // Copy the (mutable) settings objects, so that the dataset does not share them with the observation set
    const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping > dependentVariableBookkeeping =
            observationSet->getDependentVariableBookkeeping( );
    const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings = observationSet->getAncillarySettings( );
    return observationDataset.addObservationSet(
            observationSet->getObservableType( ),
            observationSet->getLinkEnds( ),
            observationSet->getObservationsReference( ),
            observationSet->getObservationTimesReference( ),
            observationSet->getReferenceLinkEnd( ),
            observationSet->getObservationsDependentVariablesReference( ),
            dependentVariableBookkeeping == nullptr ? nullptr : dependentVariableBookkeeping->clone( ),
            ancillarySettings == nullptr ? nullptr : std::make_shared< ObservationAncillarySimulationSettings >( *ancillarySettings ),
            observationSet->getWeightsReference( ),
            observationSet->getResidualsReference( ) );
}

//! Create an observation dataset from deprecated single observation sets, with one dataset set per single observation set (in order).
template< typename ObservationScalarType, typename TimeType >
std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > createObservationDataset(
        const std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >& observationSets )
{
    const std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > observationDataset =
            std::make_shared< ObservationDataset< ObservationScalarType, TimeType > >( );
    for( const auto& observationSet : observationSets )
    {
        addSingleObservationSetToDataset( *observationDataset, observationSet );
    }
    return observationDataset;
}

//! Create an observation dataset from a deprecated single observation set.
template< typename ObservationScalarType, typename TimeType >
std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > createObservationDataset(
        const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > >& observationSet )
{
    return createObservationDataset< ObservationScalarType, TimeType >(
            std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >( { observationSet } ) );
}

//! Create an observation dataset from a deprecated observation collection; dataset set ids follow the collection's set order.
template< typename ObservationScalarType, typename TimeType >
std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > createObservationDataset(
        const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >& observationCollection )
{
    if( observationCollection == nullptr )
    {
        throw std::runtime_error( "Error when creating observation dataset, input observation collection is null." );
    }
    return createObservationDataset< ObservationScalarType, TimeType >( observationCollection->getSingleObservationSets( ) );
}

//! Create a deprecated single observation set from one set of an observation dataset; only the diagonal weights are retained.
template< typename ObservationScalarType, typename TimeType >
std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > createSingleObservationSet(
        const std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > >& observationDataset,
        const unsigned int setId = 0 )
{
    if( observationDataset == nullptr )
    {
        throw std::runtime_error( "Error when creating single observation set, input observation dataset is null." );
    }
    for( const unsigned int observationId : observationDataset->getObservationIdsForSet( setId ) )
    {
        if( !observationDataset->getObservationRow( observationId ).isActive_ )
        {
            throw std::runtime_error(
                    "Error when creating single observation set from observation dataset, rejected observations are not supported by "
                    "SingleObservationSet." );
        }
    }

    const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = observationDataset->getObservationSetMetadata( setId );
    const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping =
            observationDataset->getDependentVariableBookkeeping( metadata.dependentVariableLayoutId_ );
    const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings =
            observationDataset->getAncillarySettings( metadata.ancillarySettingsId_ );
    const std::vector< Eigen::VectorXd > dependentVariables = observationDataset->getDependentVariablesForSet( setId );
    if( !dependentVariables.empty( ) && dependentVariableBookkeeping == nullptr )
    {
        throw std::runtime_error(
                "Error when creating single observation set from observation dataset, dependent variables without dependent-variable "
                "bookkeeping are not supported by SingleObservationSet." );
    }
    return std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
            metadata.observableType_,
            observationDataset->getLinkDefinition( metadata.linkDefinitionId_ ),
            observationDataset->getObservationsForSet( setId ),
            observationDataset->getObservationTimesForSet( setId ),
            metadata.referenceLinkEnd_,
            dependentVariables,
            dependentVariableBookkeeping == nullptr ? nullptr : dependentVariableBookkeeping->clone( ),
            ancillarySettings == nullptr ? nullptr : std::make_shared< ObservationAncillarySimulationSettings >( *ancillarySettings ),
            observationDataset->getWeightsForSet( setId ),
            observationDataset->getResidualsForSet( setId ) );
}

//! Create a deprecated observation collection from an observation dataset; only the diagonal weights are retained.
template< typename ObservationScalarType, typename TimeType >
std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > createObservationCollection(
        const std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > >& observationDataset )
{
    if( observationDataset == nullptr )
    {
        throw std::runtime_error( "Error when creating observation collection, input observation dataset is null." );
    }
    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > observationSets;
    for( const unsigned int setId : observationDataset->getSetIdsInObservationVectorOrder( ) )
    {
        observationSets.push_back( createSingleObservationSet( observationDataset, setId ) );
    }
    return std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >( observationSets );
}

//! Write the residuals and dependent variables of an observation dataset back into the collection it was created from.
/*!
 * The dataset must have been created from the collection with createObservationDataset (so that its set ids follow the
 * collection's set order), and its sets must still contain the same observations. Residuals and dependent variables are the
 * only per-observation data that the dataset-based computations modify.
 */
template< typename ObservationScalarType, typename TimeType >
void writeObservationDatasetToCollection(
        const ObservationDataset< ObservationScalarType, TimeType >& observationDataset,
        const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >& observationCollection )
{
    const std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > observationSets =
            observationCollection->getSingleObservationSets( );
    if( observationSets.size( ) != observationDataset.getNumberOfObservationSets( ) )
    {
        throw std::runtime_error(
                "Error when writing observation dataset to observation collection, number of observation sets is inconsistent." );
    }

    for( unsigned int setId = 0; setId < observationSets.size( ); setId++ )
    {
        if( observationSets.at( setId )->getObservableType( ) != observationDataset.getObservationSetMetadata( setId ).observableType_ ||
            observationSets.at( setId )->getObservationTimesReference( ) != observationDataset.getObservationTimesForSet( setId ) )
        {
            throw std::runtime_error(
                    "Error when writing observation dataset to observation collection, observation sets are inconsistent." );
        }

        observationSets.at( setId )->setResiduals( observationDataset.getResidualsForSet( setId ) );
        std::vector< Eigen::VectorXd > dependentVariables = observationDataset.getDependentVariablesForSet( setId );
        if( !dependentVariables.empty( ) )
        {
            observationSets.at( setId )->setObservationsDependentVariables( dependentVariables );
        }
    }
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATION_COLLECTION_CONVERSION_H

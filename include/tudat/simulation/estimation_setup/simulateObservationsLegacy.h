/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SIMULATEOBSERVATIONSLEGACY_H
#define TUDAT_SIMULATEOBSERVATIONSLEGACY_H

#include "tudat/simulation/estimation_setup/observationCollectionConversion.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"

namespace tudat
{

namespace simulation_setup
{

//! Simulate observations, returning a (deprecated) observation collection; see simulateObservationDataset.
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > simulateObservations(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > >& observationsToSimulate,
        const std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< ObservationScalarType, TimeType > > >&
                observationSimulators,
        const SystemOfBodies& bodies )
{
    return observation_models::createObservationCollection< ObservationScalarType, TimeType >(
            simulateObservationDataset< ObservationScalarType, TimeType >( observationsToSimulate, observationSimulators, bodies ) );
}

//! Create a (deprecated) observation collection from existing observations; see setExistingObservationDataset.
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > setExistingObservations(
        const std::map< observation_models::ObservableType,
                        std::pair< observation_models::LinkEnds,
                                   std::pair< std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >,
                                              std::vector< TimeType > > > > observationsInput,
        const observation_models::LinkEndType referenceLinkEnd,
        const std::map< observation_models::ObservableType, std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > >
                ancillarySettings = std::map< observation_models::ObservableType,
                                              std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > >( ) )
{
    return observation_models::createObservationCollection< ObservationScalarType, TimeType >(
            setExistingObservationDataset< ObservationScalarType, TimeType >( observationsInput, referenceLinkEnd, ancillarySettings ) );
}

//! Create observation simulation settings from a (deprecated) observation collection; see
//! getObservationSimulationSettingsFromObservationDataset.
template< typename ObservationScalarType = double, typename TimeType = double >
std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< TimeType > > >
getObservationSimulationSettingsFromObservations(
        const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observedObservationCollection,
        const SystemOfBodies& bodies )
{
    return getObservationSimulationSettingsFromObservationDataset< ObservationScalarType, TimeType >(
            observation_models::createObservationDataset( observedObservationCollection ), bodies );
}

//! Compute residuals and dependent variables for a (deprecated) observation collection, storing them in the collection; see the
//! ObservationDataset overload.
template< typename ObservationScalarType = double, typename TimeType = double >
void computeResidualsAndDependentVariables(
        const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< ObservationScalarType, TimeType > > >&
                observationSimulators,
        const SystemOfBodies& bodies )
{
    const std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > observationDataset =
            observation_models::createObservationDataset( observationCollection );
    computeResidualsAndDependentVariables< ObservationScalarType, TimeType >( observationDataset, observationSimulators, bodies );
    observation_models::writeObservationDatasetToCollection( *observationDataset, observationCollection );
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_SIMULATEOBSERVATIONSLEGACY_H

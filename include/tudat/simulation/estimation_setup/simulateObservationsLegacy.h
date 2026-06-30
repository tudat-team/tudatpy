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

#include "tudat/simulation/estimation_setup/observationCollection.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"

namespace tudat
{

namespace simulation_setup
{

template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > >
simulateObservationsWithCheckAndLinkEndIdOutput(
        const std::vector< TimeType >& observationTimes,
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel,
        const observation_models::LinkEndType referenceLinkEnd,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > > linkViabilityCalculators =
                std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > >( ),
        const std::function< Eigen::VectorXd( const double ) > noiseFunction = nullptr,
        const std::shared_ptr< ObservationDependentVariableCalculator > dependentVariableCalculator = nullptr,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr )
{
    return observation_models::createSingleObservationSet< ObservationScalarType, TimeType >(
            simulateObservationDatasetWithCheckAndLinkEndIdOutput< ObservationSize, ObservationScalarType, TimeType >(
                    observationTimes,
                    observationModel,
                    referenceLinkEnd,
                    linkViabilityCalculators,
                    noiseFunction,
                    dependentVariableCalculator,
                    ancillarySettings ) );
}

template< typename ObservationScalarType = double, typename TimeType = double, int ObservationSize = 1 >
std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > simulatePerArcSingleObservationSet(
        const std::shared_ptr< PerArcObservationSimulationSettings< TimeType > > observationsToSimulate,
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel,
        const SystemOfBodies& bodies )
{
    return observation_models::createSingleObservationSet< ObservationScalarType, TimeType >(
            simulatePerArcObservationDataset< ObservationScalarType, TimeType, ObservationSize >(
                    observationsToSimulate, observationModel, bodies ) );
}

template< typename ObservationScalarType = double, typename TimeType = double, int ObservationSize = 1 >
std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > simulateSingleObservationSet(
        const std::shared_ptr< ObservationSimulationSettings< TimeType > > observationsToSimulate,
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel,
        const SystemOfBodies& bodies )
{
    return observation_models::createSingleObservationSet< ObservationScalarType, TimeType >(
            simulateSingleObservationDataset< ObservationScalarType, TimeType, ObservationSize >(
                    observationsToSimulate, observationModel, bodies ) );
}

template< typename ObservationScalarType = double, typename TimeType = double, int ObservationSize = 1 >
std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > simulateSingleObservationSet(
        const std::shared_ptr< ObservationSimulationSettings< TimeType > > observationsToSimulate,
        const std::shared_ptr< observation_models::ObservationSimulator< ObservationSize, ObservationScalarType, TimeType > >
                observationSimulator,
        const SystemOfBodies& bodies )
{
    return observation_models::createSingleObservationSet< ObservationScalarType, TimeType >(
            simulateSingleObservationDataset< ObservationScalarType, TimeType, ObservationSize >(
                    observationsToSimulate, observationSimulator, bodies ) );
}

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

template< typename ObservationScalarType = double, typename TimeType = double >
std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< TimeType > > >
getObservationSimulationSettingsFromObservations(
        std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observedObservationCollection,
        const SystemOfBodies& bodies )
{
    return getObservationSimulationSettingsFromObservationDataset< ObservationScalarType, TimeType >(
            observedObservationCollection->getObservationDataset( ), bodies );
}

template< typename ObservationScalarType = double, typename TimeType = double >
void computeResidualsAndDependentVariables(
        std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< ObservationScalarType, TimeType > > >&
                observationSimulators,
        const SystemOfBodies& bodies )
{
    computeResidualsAndDependentVariables< ObservationScalarType, TimeType >(
            observationCollection->getObservationDataset( ), observationSimulators, bodies );
}

template< typename ObservationScalarType = double, typename TimeType = double >
void estimateTimeBiasPerSet(
        const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const Eigen::VectorXd& timePartials,
        std::vector< double >& timeBiases,
        Eigen::VectorXd& correctedResiduals )
{
    estimateTimeBiasPerSet< ObservationScalarType, TimeType >(
            observationCollection->getObservationDataset( ), timePartials, timeBiases, correctedResiduals );
}

template< typename ObservationScalarType = double, typename TimeType = double >
void estimateTimeBiasAndPolynomialFitPerSet(
        const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const Eigen::VectorXd& timePartials,
        std::vector< double >& timeBiases,
        std::vector< Eigen::VectorXd >& polynomialCoefficientsList,
        Eigen::VectorXd& correctedResiduals )
{
    estimateTimeBiasAndPolynomialFitPerSet< ObservationScalarType, TimeType >(
            observationCollection->getObservationDataset( ), timePartials, timeBiases, polynomialCoefficientsList, correctedResiduals );
}

template< typename ObservationScalarType = double, typename TimeType = double >
void getResidualStatistics(
        const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        Eigen::VectorXd& startTimes,
        Eigen::VectorXd& durations,
        Eigen::VectorXd& meanValues,
        Eigen::VectorXd& rmsValues )
{
    getResidualStatistics< ObservationScalarType, TimeType >(
            observationCollection->getObservationDataset( ), startTimes, durations, meanValues, rmsValues );
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_SIMULATEOBSERVATIONSLEGACY_H

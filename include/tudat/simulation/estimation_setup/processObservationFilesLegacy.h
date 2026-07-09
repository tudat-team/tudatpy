/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROCESSOBSERVATIONFILESLEGACY_H
#define TUDAT_PROCESSOBSERVATIONFILESLEGACY_H

#include "tudat/simulation/estimation_setup/observationCollection.h"
#include "tudat/simulation/estimation_setup/compressDopplerObservationCollection.h"
#include "tudat/simulation/estimation_setup/processOdfFile.h"
#include "tudat/simulation/estimation_setup/processPsfFile.h"
#include "tudat/simulation/estimation_setup/processTrackingTxtFile.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType = double, typename TimeType = Time >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > createOdfObservedObservationCollection(
        std::shared_ptr< ProcessedOdfFileContents< TimeType > > processedOdfFileContents,
        std::vector< observation_models::ObservableType > observableTypesToProcess = std::vector< observation_models::ObservableType >( ),
        std::pair< TimeType, TimeType > startAndEndTimesToProcess = std::make_pair< TimeType, TimeType >( TUDAT_NAN, TUDAT_NAN ),
        const bool allowDuplicateObservationsWithinObservationSet = true )
{
    return observation_models::createObservationCollection< ObservationScalarType, TimeType >(
            createOdfObservedObservationDataset< ObservationScalarType, TimeType >( processedOdfFileContents,
                                                                                    observableTypesToProcess,
                                                                                    startAndEndTimesToProcess,
                                                                                    allowDuplicateObservationsWithinObservationSet ) );
}

template< typename ObservationScalarType = double, typename TimeType = Time >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > >
createOdfObservedObservationCollectionFromFile( simulation_setup::SystemOfBodies& bodies,
                                                const std::vector< std::string >& odfFileNames,
                                                const std::string& targetName,
                                                const bool verboseOutput = true,
                                                const std::map< std::string, Eigen::Vector3d >& earthFixedGroundStationPositions =
                                                        simulation_setup::getApproximateDsnGroundStationPositions( ),
                                                const bool allowDuplicateObservationsWithinObservationSet = true )
{
    return observation_models::createObservationCollection< ObservationScalarType, TimeType >(
            createOdfObservedObservationDatasetFromFile< ObservationScalarType, TimeType >(
                    bodies,
                    odfFileNames,
                    targetName,
                    verboseOutput,
                    earthFixedGroundStationPositions,
                    allowDuplicateObservationsWithinObservationSet ) );
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::map< ObservableType, std::map< LinkEnds, std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > > >
createTrackingTxtFileObservationSets(
        const std::shared_ptr< observation_models::ProcessedTrackingTxtFileContents< ObservationScalarType, TimeType > >
                processedTrackingTxtFileContents,
        std::vector< ObservableType > observableTypesToProcess = std::vector< ObservableType >( ),
        const ObservationAncillarySimulationSettings& ancillarySettings = ObservationAncillarySimulationSettings( ) )
{
    return observation_models::createObservationCollection< ObservationScalarType, TimeType >(
                   createTrackingTxtFileObservationDataset< ObservationScalarType, TimeType >(
                           processedTrackingTxtFileContents, observableTypesToProcess, ancillarySettings ) )
            ->getObservationsSets( );
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > createTrackingTxtFilesObservationCollection(
        const std::vector< std::shared_ptr< observation_models::ProcessedTrackingTxtFileContents< ObservationScalarType, TimeType > > >
                processedTrackingTxtFileContents,
        std::vector< ObservableType > observableTypesToProcess = std::vector< ObservableType >( ),
        const ObservationAncillarySimulationSettings& ancillarySettings = ObservationAncillarySimulationSettings( ) )
{
    return observation_models::createObservationCollection< ObservationScalarType, TimeType >(
            createTrackingTxtFilesObservationDataset< ObservationScalarType, TimeType >(
                    processedTrackingTxtFileContents, observableTypesToProcess, ancillarySettings ) );
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > createTrackingTxtFileObservationCollection(
        const std::shared_ptr< observation_models::ProcessedTrackingTxtFileContents< ObservationScalarType, TimeType > >
                processedTrackingTxtFileContents,
        std::vector< ObservableType > observableTypesToProcess = std::vector< ObservableType >( ),
        const ObservationAncillarySimulationSettings& ancillarySettings = ObservationAncillarySimulationSettings( ) )
{
    return observation_models::createObservationCollection< ObservationScalarType, TimeType >(
            createTrackingTxtFileObservationDataset< ObservationScalarType, TimeType >(
                    processedTrackingTxtFileContents, observableTypesToProcess, ancillarySettings ) );
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > createTrackingTxtFileObservationCollection(
        const std::shared_ptr< input_output::TrackingTxtFileContents > rawTrackingTxtFileContents,
        const std::string spacecraftName,
        const std::vector< ObservableType > observableTypesToProcess = std::vector< ObservableType >( ),
        const std::map< std::string, Eigen::Vector3d > earthFixedGroundStationPositions =
                simulation_setup::getCombinedApproximateGroundStationPositions( ),
        const ObservationAncillarySimulationSettings& ancillarySettings = ObservationAncillarySimulationSettings( ) )
{
    auto processedTrackingTxtFileContents =
            std::make_shared< observation_models::ProcessedTrackingTxtFileContents< ObservationScalarType, TimeType > >(
                    rawTrackingTxtFileContents, spacecraftName, earthFixedGroundStationPositions );

    return createTrackingTxtFileObservationCollection( processedTrackingTxtFileContents, observableTypesToProcess, ancillarySettings );
}

template< typename ObservationScalarType = double, typename TimeType = Time >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > >
createMultiStationIfmsObservedObservationCollectionFromFiles(
        const std::vector< std::string >& ifmsFileNames,
        simulation_setup::SystemOfBodies& bodies,
        const std::string& targetName,
        const std::vector< std::string >& groundStationNames,
        const FrequencyBands& receptionBand,
        const FrequencyBands& transmissionBand,
        const bool applyTroposphereCorrection = true,
        const std::map< std::string, Eigen::Vector3d >& earthFixedGroundStationPositions =
                simulation_setup::getCombinedApproximateGroundStationPositions( ),
        const bool filterInvalidLines = true )
{
    return observation_models::createObservationCollection< ObservationScalarType, TimeType >(
            createMultiStationIfmsObservedObservationDatasetFromFiles< ObservationScalarType, TimeType >( ifmsFileNames,
                                                                                                          bodies,
                                                                                                          targetName,
                                                                                                          groundStationNames,
                                                                                                          receptionBand,
                                                                                                          transmissionBand,
                                                                                                          applyTroposphereCorrection,
                                                                                                          earthFixedGroundStationPositions,
                                                                                                          filterInvalidLines ) );
}

template< typename ObservationScalarType = double, typename TimeType = Time >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > >
createIfmsObservedObservationCollectionFromFiles( const std::vector< std::string >& ifmsFileNames,
                                                  simulation_setup::SystemOfBodies& bodies,
                                                  const std::string& targetName,
                                                  const std::string& groundStationName,
                                                  const FrequencyBands& receptionBand,
                                                  const FrequencyBands& transmissionBand,
                                                  const bool applyTroposphereCorrection = true,
                                                  const std::map< std::string, Eigen::Vector3d >& earthFixedGroundStationPositions =
                                                          simulation_setup::getCombinedApproximateGroundStationPositions( ),
                                                  const bool filterInvalidLines = true )
{
    return observation_models::createObservationCollection< ObservationScalarType, TimeType >(
            createIfmsObservedObservationDatasetFromFiles< ObservationScalarType, TimeType >( ifmsFileNames,
                                                                                              bodies,
                                                                                              targetName,
                                                                                              groundStationName,
                                                                                              receptionBand,
                                                                                              transmissionBand,
                                                                                              applyTroposphereCorrection,
                                                                                              earthFixedGroundStationPositions,
                                                                                              filterInvalidLines ) );
}

template< typename ObservationScalarType = double, typename TimeType = Time >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > >
createFdetsObservedObservationCollectionFromRawContents( const std::shared_ptr< input_output::TrackingTxtFileContents > fdetsFileContents,
                                                         const double& baseFrequency,
                                                         const std::string& targetName,
                                                         const std::string& transmittingStationName,
                                                         const std::string& receivingStationName,
                                                         const FrequencyBands& receptionBand,
                                                         const FrequencyBands& transmissionBand,
                                                         const std::map< std::string, Eigen::Vector3d >& earthFixedGroundStationPositions )
{
    return observation_models::createObservationCollection< ObservationScalarType, TimeType >(
            createFdetsObservedObservationDatasetFromRawContents< ObservationScalarType, TimeType >( fdetsFileContents,
                                                                                                     baseFrequency,
                                                                                                     targetName,
                                                                                                     transmittingStationName,
                                                                                                     receivingStationName,
                                                                                                     receptionBand,
                                                                                                     transmissionBand,
                                                                                                     earthFixedGroundStationPositions ) );
}

template< typename ObservationScalarType = double, typename TimeType = Time >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > >
createFdetsObservedObservationCollectionFromFile( const std::string& fdetsFileName,
                                                  const double& baseFrequency,
                                                  input_output::FdetDateFormat dateFormat,
                                                  const std::string& targetName,
                                                  const std::string& transmittingStationName,
                                                  const std::string& receivingStationName,
                                                  const FrequencyBands& receptionBand,
                                                  const FrequencyBands& transmissionBand,
                                                  const std::map< std::string, Eigen::Vector3d >& earthFixedGroundStationPositions =
                                                          simulation_setup::getCombinedApproximateGroundStationPositions( ) )
{
    return observation_models::createObservationCollection< ObservationScalarType, TimeType >(
            createFdetsObservedObservationDatasetFromFile< ObservationScalarType, TimeType >( fdetsFileName,
                                                                                              baseFrequency,
                                                                                              dateFormat,
                                                                                              targetName,
                                                                                              transmittingStationName,
                                                                                              receivingStationName,
                                                                                              receptionBand,
                                                                                              transmissionBand,
                                                                                              earthFixedGroundStationPositions ) );
}

template< typename ObservationScalarType = double, typename TimeType = Time >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > >
createFdetsObservedObservationCollectionFromFile( const std::string& fdetsFileName,
                                                  const double& baseFrequency,
                                                  const std::vector< std::string >& columnTypes,
                                                  const std::string& targetName,
                                                  const std::string& transmittingStationName,
                                                  const std::string& receivingStationName,
                                                  const FrequencyBands& receptionBand,
                                                  const FrequencyBands& transmissionBand,
                                                  const std::map< std::string, Eigen::Vector3d >& earthFixedGroundStationPositions =
                                                          simulation_setup::getCombinedApproximateGroundStationPositions( ) )
{
    return observation_models::createObservationCollection< ObservationScalarType, TimeType >(
            createFdetsObservedObservationDatasetFromFile< ObservationScalarType, TimeType >( fdetsFileName,
                                                                                              baseFrequency,
                                                                                              columnTypes,
                                                                                              targetName,
                                                                                              transmittingStationName,
                                                                                              receivingStationName,
                                                                                              receptionBand,
                                                                                              transmissionBand,
                                                                                              earthFixedGroundStationPositions ) );
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::map< ObservableType, std::map< LinkEnds, std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > > >
createPsfFileObservationSets( const input_output::psf::RawPsfFileContents& psfFileContents,
                              const PsfFileObservationConversionSettings& conversionSettings )
{
    return createObservationCollection< ObservationScalarType, TimeType >(
                   createPsfFileObservationDataset< ObservationScalarType, TimeType >( psfFileContents, conversionSettings ) )
            ->getObservationsSets( );
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > createPsfFileObservationCollection(
        const input_output::psf::RawPsfFileContents& psfFileContents,
        const PsfFileObservationConversionSettings& conversionSettings )
{
    return createObservationCollection< ObservationScalarType, TimeType >(
            createPsfFileObservationDataset< ObservationScalarType, TimeType >( psfFileContents, conversionSettings ) );
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > createPsfFileObservationCollection(
        const std::string& psfFile,
        const PsfFileObservationConversionSettings& conversionSettings )
{
    return createObservationCollection< ObservationScalarType, TimeType >(
            createPsfFileObservationDataset< ObservationScalarType, TimeType >( psfFile, conversionSettings ) );
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_PROCESSOBSERVATIONFILESLEGACY_H

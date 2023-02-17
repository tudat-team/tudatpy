/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PARSEODFFILE_H
#define TUDAT_PARSEODFFILE_H

#include "tudat/basics/utilities.h"
#include "tudat/io/readOdfFile.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/simulation/estimation_setup/observations.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/math/interpolators/lookupScheme.h"
#include "tudat/math/quadrature/trapezoidQuadrature.h"

namespace tudat
{

namespace orbit_determination
{

observation_models::ObservableType getObservableTypeForOdfId( const int odfId );

std::string getStationNameFromStationId ( const int networkId, const int stationId );

class ProcessedOdfFileSingleLinkData
{
public:

    ProcessedOdfFileSingleLinkData( ){ }

    virtual ~ProcessedOdfFileSingleLinkData( ){ }

    std::vector< double > observationTimes_;
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > observableValues_;
    std::vector< double > receiverDownlinkDelays_;

    std::vector< int > downlinkBandIds_;
    std::vector< int > uplinkBandIds_;
    std::vector< int > referenceBandIds_;

    std::vector< std::string > originFiles_;

    observation_models::ObservableType observableType_;

    std::string transmittingStation_;
    std::string receivingStation_;

    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > getUnprocessedObservables( )
    {
        return utilities::createMapFromVectors( observationTimes_, observableValues_ );
    }

    virtual std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > getProcessedObservables( )
    {
        return utilities::createMapFromVectors( observationTimes_, observableValues_ );
    }

    std::vector< double > getObservationTimesUtc (  )
    {
        return observationTimes_;
    }

    std::vector< double > getObservationTimesTdb (
            const simulation_setup::SystemOfBodies& bodies )
    {
        earth_orientation::TerrestrialTimeScaleConverter timeScaleConverter =
                earth_orientation::TerrestrialTimeScaleConverter( );

        std::vector< Eigen::Vector3d > earthFixedPositions;
        for ( unsigned int i = 0; i < observationTimes_.size( ); ++i )
        {
            earthFixedPositions.push_back(
                    bodies.getBody( "Earth" )->getGroundStation( receivingStation_ )->getStateInPlanetFixedFrame< double, double >(
                            observationTimes_.at( i ) ).segment( 0, 3 )
                            );
        }

        return timeScaleConverter.getCurrentTimes(
                basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, observationTimes_,
                earthFixedPositions );
    }

};


class ProcessedOdfFileDopplerData: public ProcessedOdfFileSingleLinkData
{
public:

    ~ProcessedOdfFileDopplerData( ){ }

    std::vector< int > receiverChannels_;
    std::vector< double > referenceFrequencies_;
    std::vector< double > countInterval_;
    std::vector< double > uplinkDelays_;
    std::vector< bool > receiverRampingFlags_;

    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > getProcessedObservables( )
    {
        throw std::runtime_error("Getting range rate observables not implemented");
        return utilities::createMapFromVectors( observationTimes_, observableValues_ );
    }

    std::map< double, bool > getReceiverRampingFlags( )
    {
        return utilities::createMapFromVectors( observationTimes_, receiverRampingFlags_ );
    }

    std::map< double, double > getReferenceFrequencies( )
    {
        return utilities::createMapFromVectors( observationTimes_, referenceFrequencies_ );
    }

    std::map< double, double > getCountInterval( )
    {
        return utilities::createMapFromVectors( observationTimes_, countInterval_ );
    }
};

// TODO: test computation of frequencies and integral
class RampedReferenceFrequencyInterpolator
{
public:
    RampedReferenceFrequencyInterpolator(
            std::vector< std::shared_ptr< input_output::OdfRampBlock > > rampBlock )
    {
        for( unsigned int i = 0; i < rampBlock.size( ); i++ )
        {
            startTimes_.push_back( rampBlock.at( i )->getRampStartTime( ) );
            endTimes_.push_back( rampBlock.at( i )->getRampEndTime( ) );
            rampRates_.push_back( rampBlock.at( i )->getRampRate( ) );
            startFrequencies_.push_back( rampBlock.at( i )->getRampStartFrequency( ) );
        }

        startTimeLookupScheme_ = std::make_shared<
                interpolators::HuntingAlgorithmLookupScheme< double > >(
                startTimes_ );
    }

    RampedReferenceFrequencyInterpolator(
            const std::vector< double >& startTimes_,
            const std::vector< double >& endTimes_,
            const std::vector< double >& rampRates_,
            const std::vector< double >& startFrequency_ ):
            startTimes_( startTimes_ ), endTimes_( endTimes_ ), rampRates_( rampRates_ ), startFrequencies_( startFrequency_ )
    {
        startTimeLookupScheme_ = std::make_shared<
                interpolators::HuntingAlgorithmLookupScheme< double > >(
                startTimes_ );
    }


    double getCurrentReferenceFrequencyIntegral( const double quadratureStartTime, const double quadratureEndTime )
    {
        std::vector< double > quadratureTimes;
        std::vector < double > quadratureFrequencies;

        // Point corresponding to first partial ramp
        quadratureTimes.push_back ( quadratureStartTime );
        quadratureFrequencies.push_back ( getCurrentReferenceFrequency( quadratureStartTime ) );

        // Points corresponding to full ramps
        for( unsigned int i = 1; i < startTimes_.size( ) && startTimes_.at( i ) < quadratureEndTime; i++ )
        {
            quadratureTimes.push_back( startTimes_.at( i ) );
            quadratureFrequencies.push_back( startFrequencies_.at( i ) );
        }

        // Point corresponding to final partial ramp
        quadratureTimes.push_back ( quadratureEndTime );
        quadratureFrequencies.push_back ( getCurrentReferenceFrequency( quadratureEndTime ) );

        return numerical_quadrature::performTrapezoidalQuadrature( quadratureTimes, quadratureFrequencies );
    }

    double getCurrentReferenceFrequency( const double lookupTime )
    {
        int lowerNearestNeighbour = startTimeLookupScheme_->findNearestLowerNeighbour( lookupTime );

        if( lookupTime > endTimes_.at( lowerNearestNeighbour ) || lookupTime < startTimes_.at ( lowerNearestNeighbour ) )
        {
            throw std::runtime_error(
                    "Error when interpolating ODF ramp reference frequency: look up time (" + std::to_string( lookupTime ) +
                    ") is outside the ramp table interval (" + std::to_string( startTimes_.at( 0 ) ) + " to " +
                    std::to_string( startTimes_.back( ) ) + ")." );
        }

        return startFrequencies_.at( lowerNearestNeighbour ) +
               rampRates_.at( lowerNearestNeighbour ) * ( lookupTime - startTimes_.at( lowerNearestNeighbour ) );
    }

    std::vector< double > getStartTimes ( )
    {
        return startTimes_;
    }

    std::vector< double > getEndTimes ( )
    {
        return endTimes_;
    }

    std::vector< double > getRampRates ( )
    {
        return rampRates_;
    }

    std::vector< double > getStartFrequencies ( )
    {
        return startFrequencies_;
    }

private:

    std::vector< double > startTimes_;
    std::vector< double > endTimes_;
    std::vector< double > rampRates_;
    std::vector< double > startFrequencies_;

    std::shared_ptr< interpolators::LookUpScheme< double > > startTimeLookupScheme_;

};


class ProcessedOdfFileContents
{
public:

    std::string spacecraftName_;

    std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds,
        std::shared_ptr< ProcessedOdfFileSingleLinkData > > > processedDataBlocks_;

    std::map< int, std::shared_ptr< RampedReferenceFrequencyInterpolator > > rampInterpolators_;
};


std::shared_ptr< RampedReferenceFrequencyInterpolator > mergeRampDataInterpolators(
        const std::vector< std::shared_ptr< RampedReferenceFrequencyInterpolator > >& interpolatorList );

void addOdfFileContentsToMergedContents(
        const observation_models::ObservableType observableType,
        std::shared_ptr< ProcessedOdfFileSingleLinkData > mergedOdfFileContents,
        std::shared_ptr< ProcessedOdfFileSingleLinkData > blockToAdd );

std::shared_ptr< ProcessedOdfFileContents > mergeOdfFileContents(
        const std::vector< std::shared_ptr< ProcessedOdfFileContents > > odfFileContents );

void addOdfDataBlockToProcessedData(
        const observation_models::ObservableType currentObservableType,
        const std::shared_ptr< input_output::OdfDataBlock > rawDataBlock,
        const std::shared_ptr< ProcessedOdfFileSingleLinkData > processedDataBlock );

observation_models::LinkEnds getLinkEndsFromOdfBlock (
        const std::shared_ptr< input_output::OdfDataBlock > dataBlock,
        std::string spacecraftName );

std::shared_ptr< ProcessedOdfFileContents > processOdfFileContents(
        const std::shared_ptr< input_output::OdfRawFileContents > rawOdfData,
        bool verbose = true );

} // namespace orbit_determination

} // namespace tudat

#endif // TUDAT_PARSEODFFILE_H

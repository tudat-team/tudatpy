#ifndef TUDAT_PARSEODFFILE_H
#define TUDAT_PARSEODFFILE_H

#include "tudat/basics/utilities.h"
#include "tudat/astro/orbit_determination/parseOdfFile.h"
#include "tudat/io/readOdfFile.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/math/interpolators/lookupScheme.h"

namespace tudat
{

namespace orbit_determination
{

observation_models::ObservableType getObservableTypeForOdfId(
        const int odfId );


class ProcessdOdfFileSingleLinkData
{
public:

    ProcessdOdfFileSingleLinkData( ){ }

    virtual ~ProcessdOdfFileSingleLinkData( ){ }

    std::vector< double > observationTimes;
    std::vector< double > observableValues;
    std::vector< double > receiverDownlinkDelay;

    std::vector< int > downlinkBand;
    std::vector< int > uplinkBand;
    std::vector< int > referenceBand;

    std::vector< std::string > originFile;

    observation_models::ObservableType observableType;

    std::string transmittingStation;
    std::string receivingStation;
    int transmitterNetworkId;

    std::map< double, double > getObservationData( )
    {
        return utilities::createMapFromVectors( observationTimes, observableValues );
    }
};

class ProcessdOdfFileDopplerData: public ProcessdOdfFileSingleLinkData
{
public:

    ~ProcessdOdfFileDopplerData( ){ }

    std::vector< int > receiverChannels;
    std::vector< double > referenceFrequency;
    std::vector< double > compressionTimes;
    std::vector< double > uplinkDelays;
    std::vector< double > reservedData;
    std::vector< bool > rampingFlag;

    std::map< double, bool > getRampFlags( )
    {
        return utilities::createMapFromVectors( observationTimes, rampingFlag );
    }

    std::map< double, double > getReferenceFrequencies( )
    {
        return utilities::createMapFromVectors( observationTimes, referenceFrequency );
    }

    std::map< double, double > getCompressionTimes( )
    {
        return utilities::createMapFromVectors( observationTimes, compressionTimes );
    }
};

class RampedReferenceFrequencyInterpolator
{
public:
    RampedReferenceFrequencyInterpolator(
            std::vector< std::shared_ptr< input_output::OdfRampBlock > > rampBlock )
    {
        for( unsigned int i = 0; i < rampBlock.size( ); i++ )
        {
            startTimes.push_back( rampBlock.at( i )->getRampStartTime( ) );
            endTimes.push_back( rampBlock.at( i )->getRampEndTime( ) );
            rampRates.push_back( rampBlock.at( i )->getRampRate( ) );
            startFrequency.push_back( rampBlock.at( i )->getRampStartFrequency( ) );
        }

        startTimeLookupScheme_ = std::make_shared<
                interpolators::HuntingAlgorithmLookupScheme< double > >(
                    startTimes );
    }

    RampedReferenceFrequencyInterpolator(
            const std::vector< double >& startTimes_,
            const std::vector< double >& endTimes_,
            const std::vector< double >& rampRates_,
            const std::vector< double >& startFrequency_ ):
        startTimes( startTimes_ ), endTimes( endTimes_ ), rampRates( rampRates_ ), startFrequency( startFrequency_ )
    {
        startTimeLookupScheme_ = std::make_shared<
                interpolators::HuntingAlgorithmLookupScheme< double > >(
                    startTimes );
    }


    double getCurrentReferenceFrequencyIntegral(
            const double lookupTime, const double integrationTime,
            bool& issueForInterpolation )
    {
        int lowerNearestNeighbour = startTimeLookupScheme_->findNearestLowerNeighbour(
                    lookupTime );

        if( lookupTime > endTimes.at( lowerNearestNeighbour ) )
        {
            issueForInterpolation = true;
        }

        return integrationTime * (
                    startFrequency.at( lowerNearestNeighbour ) +
                    rampRates.at( lowerNearestNeighbour ) * ( lookupTime - startTimes.at( lowerNearestNeighbour ) ) +
                    0.5 * rampRates.at( lowerNearestNeighbour ) * integrationTime );
    }

    double getCurrentReferenceFrequency(
            const double lookupTime,
            bool& issueForInterpolation )
    {
        int lowerNearestNeighbour = startTimeLookupScheme_->findNearestLowerNeighbour(
                    lookupTime );


        if( lookupTime > endTimes.at( lowerNearestNeighbour ) )
        {
            issueForInterpolation = true;
        }

        return  startFrequency.at( lowerNearestNeighbour ) +
                rampRates.at( lowerNearestNeighbour ) * ( lookupTime - startTimes.at( lowerNearestNeighbour ) );
    }


    std::vector< double > startTimes;
    std::vector< double > endTimes;
    std::vector< double > rampRates;
    std::vector< double > startFrequency;

    std::shared_ptr< interpolators::LookUpScheme< double > > startTimeLookupScheme_;

};



class ProcessedOdfFileContents
{
public:

    std::string spacecraftName;

    std::map< observation_models::ObservableType, std::map< std::pair< std::string, std::string >,
    std::shared_ptr< ProcessdOdfFileSingleLinkData > > > processedDataBlocks;

    std::map< int, std::shared_ptr< RampedReferenceFrequencyInterpolator > > rampInterpolators;
};


std::shared_ptr< RampedReferenceFrequencyInterpolator > mergeRampDataInterpolators(
        const std::vector< std::shared_ptr< RampedReferenceFrequencyInterpolator > >& interpolatorList );

void addOdfFileContentsToMergedContents(
        const observation_models::ObservableType observableType,
        std::shared_ptr< ProcessdOdfFileSingleLinkData > mergedOdfFileContents,
        std::shared_ptr< ProcessdOdfFileSingleLinkData > blockToAdd );

std::shared_ptr< ProcessedOdfFileContents > mergeOdfFileContents(
        const std::vector< std::shared_ptr< ProcessedOdfFileContents > > odfFileContents );

void addOdfDataBlockToParsedData(
        const observation_models::ObservableType currentObservableType,
        const std::shared_ptr< input_output::OdfDataBlock > rawDataBlock,
        const std::shared_ptr< ProcessdOdfFileSingleLinkData > processedDataBlock );

std::shared_ptr< ProcessedOdfFileContents > parseOdfFileContents(
        const std::shared_ptr< input_output::OdfRawFileContents > rawOdfData );

} // namespace orbit_determination

} // namespace tudat

#endif // TUDAT_PARSEODFFILE_H

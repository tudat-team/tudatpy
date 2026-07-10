/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References: 820-013, TRK-2-18 Tracking System Interfaces Orbit Data File Interface, Revision
 * E, 2008, JPL/DSN
 */

#ifndef TUDAT_PREPROCESSODFFILE_H
#define TUDAT_PREPROCESSODFFILE_H

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/basics/utilities.h"
#include "tudat/io/readOdfFile.h"
#include "tudat/io/trackingData.h"
#include "tudat/io/trackingSupplementaryData.h"
#include "tudat/math/interpolators/lookupScheme.h"

namespace tudat
{

namespace input_output
{

/*! Get the observable name associated with a given ODF data ID.
 *
 * @param odfId ODF data ID
 * @return Observable name
 */
std::string getObservableNameForOdfId( const input_output::OdfDataType observable_type_id );

/*! Get the frequency band name associated with a given ODF frequency band ID.
 *
 * @param odfId ODF frequency band ID
 * @return Frequency band name
 */
std::string getFrequencyBandNameForOdfId( const int odfId );

/*! Get the ground station name associated with a given ODF network and station ID.
 *
 * @param networkId ODF network ID
 * @param stationId ODF ground station ID
 * @return Ground station name
 */
std::string getStationNameFromStationId( const int networkId, const int stationId );

struct OdfAncillaryData {
    std::vector< int > frequencyBandIds_;
    int receptionReferenceFrequencyBandId_ = -1;
    double dopplerIntegrationTime_ = TUDAT_NAN;
    double dopplerReferenceFrequency_ = TUDAT_NAN;
    double sequentialRangeLowestRangingComponent_ = TUDAT_NAN;
    std::vector< double > linkEndsDelays_;

    bool operator==( const OdfAncillaryData& other ) const
    {
        return frequencyBandIds_ == other.frequencyBandIds_ &&
                receptionReferenceFrequencyBandId_ == other.receptionReferenceFrequencyBandId_ &&
                dopplerIntegrationTime_ == other.dopplerIntegrationTime_ &&
                dopplerReferenceFrequency_ == other.dopplerReferenceFrequency_ &&
                sequentialRangeLowestRangingComponent_ == other.sequentialRangeLowestRangingComponent_ &&
                linkEndsDelays_ == other.linkEndsDelays_;
    }
};

// Base class defining processed ODF data for a single observable and set of link ends.
template< typename TimeType = Time >
class ProcessedOdfFileSingleLinkData
{
public:
    /*!
     * Constructor.
     * @param observableType Observable type.
     * @param receivingStation Name of the receiving ground station.
     */
    ProcessedOdfFileSingleLinkData( std::string observableType, std::string receivingStation ):
        receivingStation_( receivingStation ), observableType_( observableType )
    {}

    // Destructor
    virtual ~ProcessedOdfFileSingleLinkData( ) {}

    // Observation times as seconds since EME1950 UTC
    std::vector< TimeType > unprocessedObservationTimes_;
    // Observation times as seconds since J2000 UTC
    std::vector< TimeType > processedObservationTimes_;

    // Value of the observables
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > observableValues_;
    // Vector with the downlink delay at the receiver for each observation
    std::vector< double > receiverDownlinkDelays_;

    // Vector with the downlink band ID for each observation
    std::vector< int > downlinkBandIds_;
    // Vector with the uplink band ID for each observation
    std::vector< int > uplinkBandIds_;
    // Vector with the frequency band ID of the reference frequency for each observation
    std::vector< int > referenceBandIds_;

    // Vector of files from which the data for the current observable and link ends was extracted
    std::vector< std::string > originFiles_;

    // Name of the transmitting ground station
    std::string transmittingStation_;
    // Name of the receiving ground station
    std::string receivingStation_;

    // Returns the observables mapped by the observation time
    std::map< TimeType, Eigen::Matrix< double, Eigen::Dynamic, 1 > > getObservables( )
    {
        return utilities::createMapFromVectors( processedObservationTimes_, observableValues_ );
    }

    // Returns a vector with the observables
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > getObservablesVector( )
    {
        return observableValues_;
    }

    // Returns a vector with the processed observation times
    std::vector< TimeType > getObservationTimesVector( )
    {
        return processedObservationTimes_;
    }

    std::pair< double, double > getTimeBounds( )
    {
        std::vector< double > obsTimesDouble;
        for( unsigned int k = 0; k < processedObservationTimes_.size( ); k++ )
        {
            obsTimesDouble.push_back( double( processedObservationTimes_[ k ] ) );
        }

        return std::make_pair( *std::min_element( obsTimesDouble.begin( ), obsTimesDouble.end( ) ),
                               *std::max_element( obsTimesDouble.begin( ), obsTimesDouble.end( ) ) );
    }

    // Returns the observable type
    std::string getObservableName( ) const
    {
        return observableType_;
    }

protected:
    std::string observableType_;
};

// Derived class defining Doppler (both 1- and n-way) data for a single set of link ends.
template< typename TimeType = Time >
class ProcessedOdfFileDopplerData : public ProcessedOdfFileSingleLinkData< TimeType >
{
public:
    /*!
     * Constructor.
     *
     * @param observableType Observable type.
     * @param receivingStation Name of the receiving ground station.
     * @param transmittingStation Name of the transmitting ground station.
     */
    ProcessedOdfFileDopplerData( std::string observableType, std::string receivingStation, std::string transmittingStation ):
        ProcessedOdfFileSingleLinkData< TimeType >( observableType, receivingStation ), transmittingStation_( transmittingStation )
    {}

    // Destructor
    ~ProcessedOdfFileDopplerData( ) {}

    // Name of the transmitting ground station
    std::string transmittingStation_;

    // Receiver channel per observation
    std::vector< int > receiverChannels_;
    // Reference frequency per observation
    std::vector< double > referenceFrequencies_;
    // Count interval per observation
    std::vector< double > countInterval_;
    // Uplink delay at the transmitting station per observation
    std::vector< double > transmitterUplinkDelays_;
    // Ramping flag indicating whether ramps should be used to replace receiver reference frequency
    // (if flag is false).
    std::vector< bool > receiverRampingFlags_;

    std::map< TimeType, bool > getReceiverRampingFlags( )
    {
        return utilities::createMapFromVectors( this->processedObservationTimes_, receiverRampingFlags_ );
    }

    std::map< TimeType, double > getReferenceFrequencies( )
    {
        return utilities::createMapFromVectors( this->processedObservationTimes_, referenceFrequencies_ );
    }

    std::vector< double > getReferenceFrequenciesVector( )
    {
        return referenceFrequencies_;
    }

    std::map< TimeType, double > getCountInterval( )
    {
        return utilities::createMapFromVectors( this->processedObservationTimes_, countInterval_ );
    }
};

// Derived class defining Sequential Range data for a single set of link ends.
template< typename TimeType = Time >
class ProcessedOdfFileSequentialRangeData : public ProcessedOdfFileSingleLinkData< TimeType >
{
public:
    /*!
     * Constructor.
     *
     * @param observableType Observable type.
     * @param receivingStation Name of the receiving ground station.
     * @param transmittingStation Name of the transmitting ground station.
     */
    ProcessedOdfFileSequentialRangeData( std::string observableType, std::string receivingStation, std::string transmittingStation ):
        ProcessedOdfFileSingleLinkData< TimeType >( observableType, receivingStation ), transmittingStation_( transmittingStation )
    {}

    // Destructor
    ~ProcessedOdfFileSequentialRangeData( ) {}

    // Name of the transmitting ground station
    std::string transmittingStation_;

    // Lowest ranging component per observation
    std::vector< int > lowestRangingComponent_;
    // Reference frequency per observation
    std::vector< double > referenceFrequency_;
    // Uplink ranging coder in-phase time offset from sample time tag per observation
    std::vector< int > uplinkRangingCoderInPhaseTimeOffset_;
    // Composite value per observation (highest ranging component and downlink offset)
    std::vector< double > composite2_;
    // Transmitting station uplink delay per observation
    std::vector< double > transmitterUplinkDelay_;

    // Functions to return processed data
    std::map< TimeType, int > getLowestRangingComponent( )
    {
        return utilities::createMapFromVectors( this->processedObservationTimes_, lowestRangingComponent_ );
    }

    std::map< TimeType, double > getReferenceFrequency( )
    {
        return utilities::createMapFromVectors( this->processedObservationTimes_, referenceFrequency_ );
    }

    std::map< TimeType, double > getUplinkRangingCoderInPhaseTimeOffset( )
    {
        return utilities::createMapFromVectors( this->processedObservationTimes_, uplinkRangingCoderInPhaseTimeOffset_ );
    }

    std::map< TimeType, double > getComposite( )
    {
        return utilities::createMapFromVectors( this->processedObservationTimes_, composite2_ );
    }

    std::map< TimeType, double > getTransmitterUplinkDelay( )
    {
        return utilities::createMapFromVectors( this->processedObservationTimes_, transmitterUplinkDelay_ );
    }
};

/*!
 * Compares two raw ODF data objects based on their start date. Used to sort ODF files.
 *
 * @param rawOdfData1 Raw ODF data object.
 * @param rawOdfData2 Raw ODF data object.
 * @return true if rawOdfData1 starts before rawOdfData2, false otherwise
 */
bool compareRawOdfDataByStartDate( std::shared_ptr< input_output::OdfRawFileContents > rawOdfData1,
                                   std::shared_ptr< input_output::OdfRawFileContents > rawOdfData2 );

/*!
 * Creates the link ends associated with a given ODF observation block.
 *
 * @param dataBlock ODF data block
 * @param spacecraftName Spacecraft name
 * @return Link ends
 */
data::PlainLinkDefinition getLinkEndsFromOdfBlock( const std::shared_ptr< input_output::OdfDataBlock > dataBlock,
                                                   const std::string& spacecraftName );

// Class containing processed ODF data.
template< typename TimeType = Time >
class PreProcessedOdfFileContents
{
public:
    /*!
     * Constructor for single raw ODF data object. Processes the raw ODF data.
     *
     * @param rawOdfDataVector Vector of multiple ODF data objects
     * @param spacecraftName Name of the spacecraft.
     * @param verbose Bool indicating whether to print warning regarding e.g. ignored data.
     * @param earthFixedGroundStationPositions Map with the position of each ground station in the corresponding planet's
     *      body-fixed frame. Positions are only used for converting the time between UTC and TDB,
     * therefore approximate positions are sufficient.
     */
    PreProcessedOdfFileContents( const std::shared_ptr< input_output::OdfRawFileContents > rawOdfData,
                                 const std::string& spacecraftName,
                                 bool verbose = true ):
        PreProcessedOdfFileContents( std::vector< std::shared_ptr< input_output::OdfRawFileContents > >{ rawOdfData },
                                     spacecraftName,
                                     verbose )
    {}

    /*!
     * Constructor for multiple ODF data objects. Processes the raw ODF data.
     *
     * @param rawOdfDataVector Vector of multiple ODF data objects
     * @param spacecraftName Name of the spacecraft.
     * @param verbose Bool indicating whether to print warning regarding e.g. ignored data.
     * @param earthFixedGroundStationPositions Map with the position of each ground station in the corresponding planet's
     *      body-fixed frame. Positions are only used for converting the time between UTC and TDB,
     * therefore approximate positions are sufficient.
     */
    PreProcessedOdfFileContents( std::vector< std::shared_ptr< input_output::OdfRawFileContents > > rawOdfDataVector,
                                 const std::string& spacecraftName,
                                 const bool verbose = true ):
        rawOdfData_( rawOdfDataVector ), spacecraftName_( spacecraftName ), verbose_( verbose )
    {
        // Sort ODF data files by date and check whether all the provided files apply to the same
        // spacecraft
        sortAndValidateOdfDataVector( rawOdfDataVector );

        // Extract and process raw ODF data
        extractMultipleRawOdfRampData( rawOdfDataVector );
        for( auto const& rawOdfData : rawOdfDataVector )
        {
            extractRawOdfOrbitData( rawOdfData );
        }
        printExtractionWarnings( );
        // Compute the processed observation times (i.e. TDB time from J2000)
        updateProcessedObservationTimes( );
    }

    // Get the name of the spacecraft to which the ODF data applies
    const std::string& getSpacecraftName( ) const
    {
        return spacecraftName_;
    }

    // Get the names of the ground stations included in the ODF files
    std::vector< std::string > getGroundStationsNames( )
    {
        std::vector< std::string > groundStations;

        for( auto const& [ observableType, linkDataBlocks ] : processedDataBlocks_ )
        {
            for( auto const& [ linkEnd, singleLinkDataBlock ] : linkDataBlocks )
            {
                for( auto const& [ linkEndId, linkEndType ] : linkEnd )
                {
                    // Check if linkEndId is a ground station
                    if( linkEndId[ 1 ] != "" && linkEndId[ 0 ] != spacecraftName_ )
                    {
                        if( !std::count( groundStations.begin( ), groundStations.end( ), linkEndId.getReferencePointName( ) ) )
                        {
                            groundStations.push_back( linkEndId[ 1 ] );
                        }
                    }
                }
            }
        }

        return groundStations;
    }

    // Get the observable types in the ODF files
    std::vector< std::string > getProcessedObservableNames( )
    {
        std::vector< std::string > observableNames;

        for( auto const& [ observableName, linkDataBlocks ] : processedDataBlocks_ )
        {
            observableNames.push_back( observableName );
        }

        return observableNames;
    }

    // Get pair of < start time, end time > of the data contained in the ODF files
    std::pair< double, double > getStartAndEndTime( ) const
    {
        // Reset variables
        double startTimeUtcSinceJ2000 = TUDAT_NAN;
        double endTimeUtcSinceJ2000 = TUDAT_NAN;

        // Loop over data
        for( auto const& [ observableType, linkDataBlocks ] : processedDataBlocks_ )
        {
            for( auto const& [ linkEnd, singleLinkDataBlock ] : linkDataBlocks )
            {
                // Extract the start and end times
                std::vector< TimeType > timeVectorTimeType = singleLinkDataBlock->processedObservationTimes_;

                std::vector< double > timeVector;
                for( auto const& timeTimeType : timeVectorTimeType )
                {
                    timeVector.push_back( double( timeTimeType ) );
                }

                if( timeVector.front( ) < startTimeUtcSinceJ2000 || std::isnan( startTimeUtcSinceJ2000 ) )
                {
                    startTimeUtcSinceJ2000 = timeVector.front( );
                }

                if( timeVector.back( ) > endTimeUtcSinceJ2000 || std::isnan( endTimeUtcSinceJ2000 ) )
                {
                    endTimeUtcSinceJ2000 = timeVector.back( );
                }
            }
        }

        return std::make_pair( startTimeUtcSinceJ2000, endTimeUtcSinceJ2000 );
    }

    // Return ODF observable types IDs (as per TRK-2-18) that were not included in the processed
    // data
    std::vector< input_output::OdfDataType > getIgnoredRawOdfObservableTypes( ) const
    {
        return ignoredRawOdfObservableTypes_;
    }

    // Return ground stations for which observations were not included in the processed data (due to
    // absence of ramp tables)
    std::vector< std::string > getIgnoredGroundStations( ) const
    {
        return ignoredGroundStations_;
    }

    // Return the ramp tables per ground station
    const std::map< std::string, std::shared_ptr< data::RampedFrequencySupplementaryData > >& getRampTables( ) const
    {
        return rampTables_;
    }

    // Return processed data
    const std::map< std::string, std::map< data::PlainLinkDefinition, std::shared_ptr< ProcessedOdfFileSingleLinkData< TimeType > > > >&
    getProcessedDataBlocks( ) const
    {
        return processedDataBlocks_;
    }

    // Return the raw ODF data
    std::vector< std::shared_ptr< input_output::OdfRawFileContents > > getRawOdfData( )
    {
        return rawOdfData_;
    }

    void defineSpacecraftAntennaId( const std::string& spacecraft, const std::string& antennaName )
    {
        std::map< std::string, std::map< data::PlainLinkDefinition, std::shared_ptr< ProcessedOdfFileSingleLinkData< TimeType > > > >
                reprocessedDataBlocks;

        for( auto const& [ observableType, linkDataBlocks ] : processedDataBlocks_ )
        {
            for( auto const& [ linkEnd, singleLinkDataBlock ] : linkDataBlocks )
            {
                data::PlainLinkDefinition newLinkEnds;
                for( auto const& [ linkEndId, linkEndType ] : linkEnd )
                {
                    if( linkEndId[ 0 ] == spacecraft )
                    {
                        newLinkEnds.push_back( std::make_pair( std::make_pair( spacecraft, antennaName ), linkEndType ) );
                    }
                    else
                    {
                        newLinkEnds.push_back( std::make_pair( linkEndId, linkEndType ) );
                    }
                }
                reprocessedDataBlocks[ observableType ][ newLinkEnds ] = singleLinkDataBlock;
            }
        }
        processedDataBlocks_ = reprocessedDataBlocks;
    }

    void defineSpacecraftAntennaId( const std::string& spacecraft,
                                    const std::string& antennaName,
                                    const std::map< double, double >& timeIntervals )
    {
        // Create lookup scheme to find closest time interval
        std::vector< double > timeIntervalStarts = utilities::createVectorFromMapKeys( timeIntervals );
        auto lookUpScheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >( timeIntervalStarts );

        // Define new data blocks
        std::map< std::string, std::map< data::PlainLinkDefinition, std::shared_ptr< ProcessedOdfFileSingleLinkData< TimeType > > > >
                reprocessedDataBlocks;

        // Iterate over all observable types
        for( auto const& [ observableType, linkDataBlocks ] : processedDataBlocks_ )
        {
            // Iterate over all link ends
            for( auto const& [ linkEnd, singleLinkDataBlock ] : linkDataBlocks )
            {
                // Get start and end times of current block
                std::pair< double, double > timeBounds = singleLinkDataBlock->getTimeBounds( );

                // Get start and end time of nearest interval
                double timeIntervalStart = TUDAT_NAN;
                if( timeBounds.first > timeIntervalStarts.at( 0 ) )
                {
                    int timeIntervalIndex = lookUpScheme_->findNearestLowerNeighbour( timeBounds.first );
                    timeIntervalStart = timeIntervalStarts.at( timeIntervalIndex );
                }
                else
                {
                    timeIntervalStart = timeIntervalStarts.at( 0 );
                }
                double timeIntervalEnd = timeIntervalStart + timeIntervals.at( timeIntervalStart );

                // Copy entire block
                if( timeIntervalStart < timeBounds.first && timeIntervalEnd > timeBounds.second )
                {
                    data::PlainLinkDefinition newLinkEnds;
                    for( auto const& [ linkEndId, linkEndType ] : linkEnd )
                    {
                        if( linkEndId[ 0 ] == spacecraft )
                        {
                            newLinkEnds.push_back( std::make_pair( std::make_pair( spacecraft, antennaName ), linkEndType ) );
                        }
                        else
                        {
                            newLinkEnds.push_back( std::make_pair( linkEndId, linkEndType ) );
                        }
                    }
                    reprocessedDataBlocks[ observableType ][ newLinkEnds ] = singleLinkDataBlock;
                }
                else if( timeIntervalEnd > timeBounds.second )
                {
                }
            }
        }
        processedDataBlocks_ = reprocessedDataBlocks;
    }

private:
    /*!
     * Checks whether the vector of ODF data is valid (i.e. all objects apply to the same
     * spacecraft), and if so, sorts the vector by the date of the ODF objets.
     *
     * @param rawOdfDataVector Vector of raw ODF objects.
     */
    void sortAndValidateOdfDataVector( std::vector< std::shared_ptr< input_output::OdfRawFileContents > >& rawOdfDataVector )
    {
        unsigned int spacecraftId = rawOdfDataVector.front( )->spacecraftId_;

        for( auto const& rawOdfData : rawOdfDataVector )
        {
            // Check if spacecraft ID is valid
            if( rawOdfData->spacecraftId_ != spacecraftId )
            {
                throw std::runtime_error(
                        "Error when creating processed ODF object from raw data: multiple "
                        "spacecraft IDs"
                        "found (" +
                        std::to_string( spacecraftId ) + " and " + std::to_string( rawOdfData->spacecraftId_ ) + ")." );
            }
        }

        std::stable_sort( rawOdfDataVector.begin( ), rawOdfDataVector.end( ), &compareRawOdfDataByStartDate );
    }

    /*!
     * Checks whether a given observation is valid. Checks if the observation time is covered by the
     * available ramp tables, for the relevant ground station(s).
     *
     * @param rawDataBlock Raw ODF data block
     * @param linkEnds Link ends to which the ODF block applies
     * @param currentObservableType Observable type
     * @return Bool indicating whether observation is valid or not
     */
    bool isObservationValid( std::shared_ptr< input_output::OdfDataBlock > rawDataBlock, const std::string& currentObservableType )
    {
        input_output::OdfDataType currentObservableId = rawDataBlock->getObservableSpecificDataBlock( )->dataType_;

        if( requiresTransmittingStation( observation_models::getObservableType( currentObservableType ) ) )
        {
            std::string transmittingStation =
                    getStationNameFromStationId( rawDataBlock->getCommonDataBlock( )->transmittingStationNetworkId_,
                                                 rawDataBlock->getCommonDataBlock( )->transmittingStationId_ );

            // Check if transmitting station is in ramp tables
            if( rampTables_.count( transmittingStation ) == 0 )
            {
                if( std::count( ignoredGroundStations_.begin( ), ignoredGroundStations_.end( ), transmittingStation ) == 0 )
                {
                    ignoredGroundStations_.push_back( transmittingStation );
                    if( verbose_ )
                    {
                        std::cerr << "Warning: ground station " << transmittingStation << " not available in ramp tables,"
                                  << " ignoring corresponding data." << std::endl;
                    }
                }
                ignoredOdfRawDataBlocks_.push_back( rawDataBlock );
                return false;
            }

            // Check if observation time is covered by ramp tables
            if( rawDataBlock->getCommonDataBlock( )->getObservableTime( ) <
                unprocessedRampStartTimesPerStation_[ transmittingStation ].front( ) )
            {
                if( verbose_ )
                {
                    noRampDataItems_[ static_cast< int >( currentObservableId ) ][ transmittingStation ].push_back(
                            rawDataBlock->getCommonDataBlock( )->getObservableTime( ) );
                }
                ignoredOdfRawDataBlocks_.push_back( rawDataBlock );
                return false;
            }
        }
        if( requiresFirstReceivingStation( observation_models::getObservableType( currentObservableType ) ) )
        {
            std::string receivingStation = getStationNameFromStationId( 0, rawDataBlock->getCommonDataBlock( )->receivingStationId_ );

            // Check if receiving station is in ramp tables
            if( rampTables_.count( receivingStation ) == 0 )
            {
                if( std::count( ignoredGroundStations_.begin( ), ignoredGroundStations_.end( ), receivingStation ) == 0 )
                {
                    ignoredGroundStations_.push_back( receivingStation );
                    if( verbose_ )
                    {
                        std::cerr << "Warning: observation of ODF type " << static_cast< int >( currentObservableId )
                                  << " not covered by ramp table of station " << receivingStation << ", ignoring it." << std::endl;
                    }
                }
                ignoredOdfRawDataBlocks_.push_back( rawDataBlock );
                return false;
            }

            // Check if observation time is covered by ramp tables
            if( rawDataBlock->getCommonDataBlock( )->getObservableTime( ) <
                        unprocessedRampStartTimesPerStation_[ receivingStation ].front( ) ||
                rawDataBlock->getCommonDataBlock( )->getObservableTime( ) >
                        unprocessedRampStartTimesPerStation_[ receivingStation ].back( ) )
            {
                if( verbose_ )
                {
                    std::cerr << "Warning: observation of ODF type " << static_cast< int >( currentObservableId )
                              << " not covered by ramp tables," << " ignoring it." << std::endl;
                }
                ignoredOdfRawDataBlocks_.push_back( rawDataBlock );
                return false;
            }
        }

        return true;
    }

    /*!
     * Extracts data from a raw ODF file, splitting it based on observable type and link ends.
     *
     * @param rawOdfData Raw ODF data object.
     */
    void extractRawOdfOrbitData( std::shared_ptr< input_output::OdfRawFileContents > rawOdfData )
    {
        // Retrieve data blocks from ODF file raw contents
        std::vector< std::shared_ptr< input_output::OdfDataBlock > > rawDataBlocks = rawOdfData->getDataBlocks( );

        // Iterate over all block of ODF file.
        for( auto const& rawDataBlock : rawDataBlocks )
        {
            // Retrieve observable type and link ends
            input_output::OdfDataType currentObservableId = rawDataBlock->getObservableSpecificDataBlock( )->dataType_;

            // Get current observable type and throw warning if not implemented
            std::string currentObservableType;
            try
            {
                currentObservableType = getObservableNameForOdfId( currentObservableId );
            }
            catch( const std::runtime_error& )
            {
                if( std::find( ignoredRawOdfObservableTypes_.begin( ), ignoredRawOdfObservableTypes_.end( ), currentObservableId ) ==
                    ignoredRawOdfObservableTypes_.end( ) )
                {
                    ignoredRawOdfObservableTypes_.push_back( currentObservableId );
                    if( verbose_ )
                    {
                        std::cerr << "Warning: processing of ODF data type " << static_cast< int >( currentObservableId )
                                  << " is not implemented, ignoring the corresponding data." << std::endl;
                    }
                }
                ignoredOdfRawDataBlocks_.push_back( rawDataBlock );
                continue;
            }

            auto linkEnds = getLinkEndsFromOdfBlock( rawDataBlock, spacecraftName_ );

            // Check if observation is valid and should be processed
            if( isObservationValid( rawDataBlock, currentObservableType ) )
            {
                // Check if data object already exists for current observable/link ends
                bool createNewObject = false;
                if( processedDataBlocks_.count( currentObservableType ) == 0 )
                {
                    createNewObject = true;
                }
                else if( processedDataBlocks_.at( currentObservableType ).count( linkEnds ) == 0 )
                {
                    createNewObject = true;
                }

                // Create new data object, if required
                if( createNewObject )
                {
                    if( currentObservableType == "DsnNWayAveragedDoppler" )
                    {
                        const std::string receivingStation =
                                getStationNameFromStationId( 0, rawDataBlock->getCommonDataBlock( )->receivingStationId_ );
                        const std::string transmittingStation =
                                getStationNameFromStationId( rawDataBlock->getCommonDataBlock( )->transmittingStationNetworkId_,
                                                             rawDataBlock->getCommonDataBlock( )->transmittingStationId_ );
                        processedDataBlocks_[ currentObservableType ][ linkEnds ] =
                                std::make_shared< ProcessedOdfFileDopplerData< TimeType > >(
                                        currentObservableType, receivingStation, transmittingStation );
                    }
                    else if( currentObservableType == "DsnNWayRange" )
                    {
                        const std::string receivingStation =
                                getStationNameFromStationId( 0, rawDataBlock->getCommonDataBlock( )->receivingStationId_ );
                        const std::string transmittingStation =
                                getStationNameFromStationId( rawDataBlock->getCommonDataBlock( )->transmittingStationNetworkId_,
                                                             rawDataBlock->getCommonDataBlock( )->transmittingStationId_ );
                        processedDataBlocks_[ currentObservableType ][ linkEnds ] =
                                std::make_shared< ProcessedOdfFileSequentialRangeData< TimeType > >(
                                        currentObservableType, receivingStation, transmittingStation );
                    }
                    else
                    {
                        std::cerr << "Warning: Observable type " << ( currentObservableType ) << " is not implemented." << std::endl;
                    }
                }
            }

            addOdfRawDataBlockToProcessedData(
                    rawDataBlock, processedDataBlocks_[ currentObservableType ][ linkEnds ], rawOdfData->fileName_ );
        }
    }

    void printExtractionWarnings( )
    {
        for( auto it : noRampDataItems_ )
        {
            for( auto it2 : it.second )
            {
                std::cerr << "Warning: observation of ODF type " << it.first << ", " << it2.second.size( )
                          << " observations with transmitting station " << it2.first
                          << " not covered by ramp table of station. These observations are ignored and not processed further."
                          << std::endl;
            }
        }
    }

    /*!
     * Add an unprocessed ODF data block to the processed data object associated with the relevant
     * observable type and link ends.
     *
     * @param rawDataBlock Raw ODF data block.
     * @param singleLinkProcessedData Processed data for the link ends to which the raw data block applies.
     * @param rawDataFileName Name of the file from which the raw data was extracted.
     */
    void addOdfRawDataBlockToProcessedData( const std::shared_ptr< input_output::OdfDataBlock > rawDataBlock,
                                            const std::shared_ptr< ProcessedOdfFileSingleLinkData< TimeType > > singleLinkProcessedData,
                                            const std::string& rawDataFileName )
    {
        // Add properties to data block if data is valid
        if( rawDataBlock->getCommonDataBlock( )->validity_ == 0 )
        {
            // Add common properties to data object
            singleLinkProcessedData->downlinkBandIds_.push_back( rawDataBlock->getCommonDataBlock( )->downlinkBandId_ );
            singleLinkProcessedData->uplinkBandIds_.push_back( rawDataBlock->getCommonDataBlock( )->uplinkBandId_ );
            singleLinkProcessedData->referenceBandIds_.push_back( rawDataBlock->getCommonDataBlock( )->referenceBandId_ );
            singleLinkProcessedData->unprocessedObservationTimes_.push_back( rawDataBlock->getCommonDataBlock( )->getObservableTime( ) );
            singleLinkProcessedData->receiverDownlinkDelays_.push_back(
                    rawDataBlock->getCommonDataBlock( )->getReceivingStationDownlinkDelay( ) );
            singleLinkProcessedData->originFiles_.push_back( rawDataFileName );

            // Add properties to data object for Doppler data
            if( singleLinkProcessedData->getObservableName( ) == "DsnNWayAveragedDoppler" )
            {
                std::shared_ptr< input_output::OdfDopplerDataBlock > odfDopplerDataBlock =
                        std::dynamic_pointer_cast< input_output::OdfDopplerDataBlock >( rawDataBlock->getObservableSpecificDataBlock( ) );
                std::shared_ptr< ProcessedOdfFileDopplerData< TimeType > > odfParsedDopplerDataBlock =
                        std::dynamic_pointer_cast< ProcessedOdfFileDopplerData< TimeType > >( singleLinkProcessedData );

                singleLinkProcessedData->observableValues_.push_back(
                        ( Eigen::Matrix< double, 1, 1 >( ) << rawDataBlock->getCommonDataBlock( )->getObservableValue( ) ).finished( ) );

                odfParsedDopplerDataBlock->countInterval_.push_back( odfDopplerDataBlock->getCompressionTime( ) );
                odfParsedDopplerDataBlock->receiverChannels_.push_back( odfDopplerDataBlock->getReceiverChannel( ) );
                odfParsedDopplerDataBlock->receiverRampingFlags_.push_back( odfDopplerDataBlock->getReceiverExciterFlag( ) );
                odfParsedDopplerDataBlock->referenceFrequencies_.push_back( odfDopplerDataBlock->getReferenceFrequency( ) );
                odfParsedDopplerDataBlock->transmitterUplinkDelays_.push_back( odfDopplerDataBlock->getTransmittingStationUplinkDelay( ) );
            }
            else if( singleLinkProcessedData->getObservableName( ) == "DsnNWayRange" )
            {
                std::shared_ptr< input_output::OdfSequentialRangeDataBlock > odfSequentialRangeDataBlock =
                        std::dynamic_pointer_cast< input_output::OdfSequentialRangeDataBlock >(
                                rawDataBlock->getObservableSpecificDataBlock( ) );
                std::shared_ptr< ProcessedOdfFileSequentialRangeData< TimeType > > odfParsedSequentialRangeDataBlock =
                        std::dynamic_pointer_cast< ProcessedOdfFileSequentialRangeData< TimeType > >( singleLinkProcessedData );

                singleLinkProcessedData->observableValues_.push_back(
                        ( Eigen::Matrix< double, 1, 1 >( ) << rawDataBlock->getCommonDataBlock( )->getObservableValue( ) ).finished( ) );

                odfParsedSequentialRangeDataBlock->composite2_.push_back( odfSequentialRangeDataBlock->getCompositeTwo( ) );
                odfParsedSequentialRangeDataBlock->lowestRangingComponent_.push_back(
                        odfSequentialRangeDataBlock->getLowestRangingComponent( ) );
                odfParsedSequentialRangeDataBlock->referenceFrequency_.push_back( odfSequentialRangeDataBlock->getReferenceFrequency( ) );
                odfParsedSequentialRangeDataBlock->uplinkRangingCoderInPhaseTimeOffset_.push_back(
                        odfSequentialRangeDataBlock->getUplinkCoderInPhaseTimeOffset( ) );
                odfParsedSequentialRangeDataBlock->transmitterUplinkDelay_.push_back(
                        odfSequentialRangeDataBlock->getTransmittingStationUplinkDelay( ) );
            }
        }
    }

    /*!
     * Extracts and merges the ramp data from the provided ODF files, creating one frequency
     * interpolator object per ground station.
     *
     * @param rawOdfDataVector Vector of raw ODF data objects.
     */
    void extractMultipleRawOdfRampData( std::vector< std::shared_ptr< input_output::OdfRawFileContents > > rawOdfDataVector )
    {
        std::map< std::string, std::vector< double > > rampRatesPerStation, startFrequenciesPerStation;

        for( auto const& rawOdfData : rawOdfDataVector )
        {
            std::map< int, std::vector< std::shared_ptr< input_output::OdfRampBlock > > > rampBlocksPerStation =
                    rawOdfData->getRampBlocks( );
            for( auto const& [ stationId, rampBlocks ] : rampBlocksPerStation )
            {
                std::string stationName = getStationNameFromStationId( 0, stationId );

                for( unsigned int j = 0; j < rampBlocks.size( ); j++ )
                {
                    // Check if zero time ramp
                    if( rampBlocks.at( j )->getRampStartTime( ) == rampBlocks.at( j )->getRampEndTime( ) )
                    {
                        continue;
                    }

                    auto& rampStartTimes = unprocessedRampStartTimesPerStation_[ stationName ];
                    auto& rampEndTimes = unprocessedRampEndTimesPerStation_[ stationName ];
                    auto& rampRates = rampRatesPerStation[ stationName ];
                    auto& startFrequencies = startFrequenciesPerStation[ stationName ];

                    // Add a connection interval only when consecutive ODF files leave a real gap.
                    if( j == 0 && !rampStartTimes.empty( ) )
                    {
                        Time previousRampEndTime = rampEndTimes.back( );
                        Time currentRampStartTime = rampBlocks.at( j )->getRampStartTime( );
                        if( previousRampEndTime < currentRampStartTime )
                        {
                            rampStartTimes.push_back( previousRampEndTime );
                            rampEndTimes.push_back( currentRampStartTime );
                            rampRates.push_back( TUDAT_NAN );
                            startFrequencies.push_back( TUDAT_NAN );
                        }
                    }

                    rampStartTimes.push_back( rampBlocks.at( j )->getRampStartTime( ) );
                    rampEndTimes.push_back( rampBlocks.at( j )->getRampEndTime( ) );
                    rampRates.push_back( rampBlocks.at( j )->getRampRate( ) );
                    startFrequencies.push_back( rampBlocks.at( j )->getRampStartFrequency( ) );
                }
            }
        }

        for( auto it = unprocessedRampStartTimesPerStation_.begin( ); it != unprocessedRampStartTimesPerStation_.end( ); ++it )
        {
            std::string stationName = it->first;

            std::vector< Time > rampStartTimesPerStationUtc =
                    computeObservationTimesUtcFromJ2000< Time >( unprocessedRampStartTimesPerStation_[ stationName ] );
            std::vector< Time > rampEndTimesPerStationUtc =
                    computeObservationTimesUtcFromJ2000< Time >( unprocessedRampEndTimesPerStation_[ stationName ] );

            rampTables_[ stationName ] = std::make_shared< data::RampedFrequencySupplementaryData >( );

            for( unsigned int i = 0; i < rampStartTimesPerStationUtc.size( ); ++i )
            {
                rampTables_[ stationName ]->addFrequencyRamp( rampStartTimesPerStationUtc.at( i ),
                                                              rampEndTimesPerStationUtc.at( i ),
                                                              rampRatesPerStation[ stationName ].at( i ),
                                                              startFrequenciesPerStation[ stationName ].at( i ) );
            }
        }
    }

    /*!
     * Goes over all the extracted ibservations and converts the observation times to TDB from
     * J2000.
     */
    void updateProcessedObservationTimes( )
    {
        // Loop over saved data and convert time to TDB wrt J2000
        for( auto const& [ observableType, linkDataBlocks ] : processedDataBlocks_ )
        {
            for( auto const& [ linkEnd, singleLinkDataBlock ] : linkDataBlocks )
            {
                singleLinkDataBlock->processedObservationTimes_ =
                        computeObservationTimesUtcFromJ2000< Time >( singleLinkDataBlock->unprocessedObservationTimes_ );
            }
        }
    }

    /*!
     * Converts UTC times from EME1950 to UTC times from J2000.
     *
     * @param observationTimesUtcFromEME1950 UTC times from EME1950.
     * @return UTC times from J2000
     */
    template< typename InputTimeType >
    std::vector< Time > computeObservationTimesUtcFromJ2000( const std::vector< InputTimeType >& observationTimesUtcFromEME1950 )
    {
        std::vector< Time > observationTimesUtcFromJ2000;

        auto EME1950ToJ2000Offset = Time( basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch< double >(
                                                  1950, 1, 1, 0, 0, 0, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
                                          physical_constants::JULIAN_DAY );

        for( unsigned int i = 0; i < observationTimesUtcFromEME1950.size( ); ++i )
        {
            observationTimesUtcFromJ2000.push_back( static_cast< Time >( observationTimesUtcFromEME1950.at( i ) ) + EME1950ToJ2000Offset );
        }

        return observationTimesUtcFromJ2000;
    }

    // Vector of raw ODF data
    std::vector< std::shared_ptr< input_output::OdfRawFileContents > > rawOdfData_;

    // Name of the spacecraft
    const std::string spacecraftName_;

    //    const std::string antennaName_;

    // Processed data mapped by observable type and link ends
    std::map< std::string, std::map< data::PlainLinkDefinition, std::shared_ptr< ProcessedOdfFileSingleLinkData< TimeType > > > >
            processedDataBlocks_;

    // Transmitting frequency objects mapped by the ground station names
    std::map< std::string, std::shared_ptr< data::RampedFrequencySupplementaryData > > rampTables_;

    // Unprocessed ramp start and end times. Used for pre-processing observations
    std::map< std::string, std::vector< Time > > unprocessedRampStartTimesPerStation_;
    std::map< std::string, std::vector< Time > > unprocessedRampEndTimesPerStation_;

    // Ignored data: either because processing of the data type has not been implemented or because
    // the required data to simulate the corresponding observable is not available in the ODF files
    // (e.g. ramp tables located in some other ODF file that was not provided). Vector keeping the
    // invalid observable types that were found in the raw ODF data
    std::vector< input_output::OdfDataType > ignoredRawOdfObservableTypes_;
    // Vector keeping the invalid ground stations that were found in the raw ODF data
    std::vector< std::string > ignoredGroundStations_;
    // Vector keeping the invalid data blocks that were found in the raw ODF data
    std::vector< std::shared_ptr< input_output::OdfDataBlock > > ignoredOdfRawDataBlocks_;

    // Flag indicating whether to print warnings
    bool verbose_;

    std::map< int, std::map< std::string, std::vector< Time > > > noRampDataItems_;

    // TODO: friend class used in unit test. Remove after processing of ODF data type 11 (1-way
    // Doppler) is implemented
    friend class PreProcessedOdfFileContentsPrivateFunctionTest;
};

template< typename TimeType = Time >
inline std::shared_ptr< PreProcessedOdfFileContents< TimeType > > processOdfData( const std::vector< std::string >& odfFileNames,
                                                                                  const std::string& spacecraftName,
                                                                                  const bool verbose = true )
{
    std::vector< std::shared_ptr< input_output::OdfRawFileContents > > rawOdfDataVector;
    for( std::string odfFile : odfFileNames )
    {
        rawOdfDataVector.push_back( std::make_shared< input_output::OdfRawFileContents >( odfFile ) );
    }

    return std::make_shared< PreProcessedOdfFileContents< TimeType > >( rawOdfDataVector, spacecraftName, verbose );
}

template< typename TimeType = Time >
inline std::shared_ptr< PreProcessedOdfFileContents< TimeType > > processOdfData( const std::string& odfFileName,
                                                                                  const std::string& spacecraftName,
                                                                                  const bool verbose = true )
{
    return processOdfData< TimeType >( std::vector< std::string >{ odfFileName }, spacecraftName, verbose );
}

template< typename TimeType = Time >
inline std::shared_ptr< PreProcessedOdfFileContents< TimeType > > processOdfData(
        const std::vector< std::shared_ptr< input_output::OdfRawFileContents > >& odfFiles,
        const std::string& spacecraftName,
        const bool verbose = true )
{
    return std::make_shared< PreProcessedOdfFileContents< TimeType > >( odfFiles, spacecraftName, verbose );
}

template< typename TimeType = Time >
inline std::shared_ptr< PreProcessedOdfFileContents< TimeType > > processOdfData(
        const std::shared_ptr< input_output::OdfRawFileContents > odfFile,
        const std::string& spacecraftName,
        const bool verbose = true )
{
    return processOdfData< TimeType >(
            std::vector< std::shared_ptr< input_output::OdfRawFileContents > >{ odfFile }, spacecraftName, verbose );
}

/*!
 * Creates the plain ODF metadata for the observations indexed by dataIndex in the provided
 * processed ODF data.
 *
 * @param odfDataContents Processed ODF data.
 * @param dataIndex Index of the observation for which to create the metadata.
 * @return ODF observation metadata.
 */
template< typename TimeType = double >
OdfAncillaryData createOdfAncillaryData( std::shared_ptr< ProcessedOdfFileSingleLinkData< TimeType > > odfDataContents,
                                         unsigned int dataIndex )
{
    if( dataIndex >= odfDataContents->unprocessedObservationTimes_.size( ) )
    {
        throw std::runtime_error(
                "Error when creating ODF data metadata: specified data index is larger "
                "than data size." );
    }

    OdfAncillaryData ancillaryData;

    const std::string currentObservableType = odfDataContents->getObservableName( );

    ancillaryData.frequencyBandIds_ =
            std::vector< int >{ odfDataContents->uplinkBandIds_.at( dataIndex ), odfDataContents->downlinkBandIds_.at( dataIndex ) };
    ancillaryData.receptionReferenceFrequencyBandId_ = odfDataContents->referenceBandIds_.at( dataIndex );

    if( auto dopplerDataBlock = std::dynamic_pointer_cast< ProcessedOdfFileDopplerData< TimeType > >( odfDataContents ) )
    {
        ancillaryData.dopplerIntegrationTime_ = dopplerDataBlock->countInterval_.at( dataIndex );
        ancillaryData.dopplerReferenceFrequency_ = dopplerDataBlock->referenceFrequencies_.at( dataIndex );

        if( currentObservableType == "DsnNWayAveragedDoppler" )
        {
            ancillaryData.linkEndsDelays_ = std::vector< double >{ dopplerDataBlock->transmitterUplinkDelays_.at( dataIndex ),
                                                                   0.0,
                                                                   dopplerDataBlock->receiverDownlinkDelays_.at( dataIndex ) };
        }
        else
        {
            ancillaryData.linkEndsDelays_ = std::vector< double >{ dopplerDataBlock->transmitterUplinkDelays_.at( dataIndex ),
                                                                   dopplerDataBlock->receiverDownlinkDelays_.at( dataIndex ) };
        }
    }
    else if( auto sequentialRangeDataBlock =
                     std::dynamic_pointer_cast< ProcessedOdfFileSequentialRangeData< TimeType > >( odfDataContents ) )
    {
        ancillaryData.sequentialRangeLowestRangingComponent_ = sequentialRangeDataBlock->lowestRangingComponent_.at( dataIndex );
        ancillaryData.linkEndsDelays_ = std::vector< double >{ sequentialRangeDataBlock->transmitterUplinkDelay_.at( dataIndex ),
                                                               0.0,
                                                               sequentialRangeDataBlock->receiverDownlinkDelays_.at( dataIndex ) };
    }
    else
    {
        throw std::runtime_error( "Error when casting ODF processed data: data type not identified." );
    }

    return ancillaryData;
}

/*!
 * Given the processed ODF data for a single set of link ends, extracts for each observation the
 * observation time, observable value, and ODF ancillary data.
 *
 * @param odfSingleLinkData Observations data for a single set of link ends (input)
 * @param observationTimes Vector of observation times (output)
 * @param observables Vector of observables (output)
 * @param ancillaryData Vector of ODF ancillary data (output)
 */
template< typename ObservationScalarType = double, typename TimeType = double >
void separateSingleLinkOdfData( std::shared_ptr< ProcessedOdfFileSingleLinkData< TimeType > > odfSingleLinkData,
                                std::vector< std::vector< TimeType > >& observationTimes,
                                std::vector< std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > >& observables,
                                std::vector< OdfAncillaryData >& ancillaryData )
{
    // Initialize vectors
    observationTimes.clear( );
    observables.clear( );
    ancillaryData.clear( );

    // Get time and observables vectors
    std::vector< TimeType > observationTimesTdb = odfSingleLinkData->getObservationTimesVector( );
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > observablesVector = odfSingleLinkData->getObservablesVector( );

    for( unsigned int i = 0; i < odfSingleLinkData->unprocessedObservationTimes_.size( ); ++i )
    {
        OdfAncillaryData currentAncillaryData = createOdfAncillaryData< TimeType >( odfSingleLinkData, i );

        bool newAncillaryData = true;

        for( unsigned int j = 0; j < ancillaryData.size( ); ++j )
        {
            if( ancillaryData.at( j ) == currentAncillaryData )
            {
                newAncillaryData = false;
                observationTimes.at( j ).push_back( observationTimesTdb.at( i ) );
                observables.at( j ).push_back( observablesVector.at( i ).template cast< ObservationScalarType >( ) );
                break;
            }
        }

        if( newAncillaryData )
        {
            observationTimes.push_back( std::vector< TimeType >{ observationTimesTdb.at( i ) } );
            observables.push_back( std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >{
                    observablesVector.at( i ).template cast< ObservationScalarType >( ) } );
            ancillaryData.push_back( currentAncillaryData );
        }
    }
}

template< typename ObservationScalarType = double, typename TimeType = Time >
void setOdfMetadataInTrackingData( std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > trackingData,
                                   const OdfAncillaryData& ancillaryData )
{
    trackingData->addAncillarySettings(
            "frequency bands",
            std::vector< std::string >{ getFrequencyBandNameForOdfId( ancillaryData.frequencyBandIds_.at( 0 ) ),
                                        getFrequencyBandNameForOdfId( ancillaryData.frequencyBandIds_.at( 1 ) ) } );
    trackingData->addAncillarySettings( "DSN reference frequency band at reception",
                                        getFrequencyBandNameForOdfId( ancillaryData.receptionReferenceFrequencyBandId_ ) );

    if( ( std::isnan( ancillaryData.dopplerIntegrationTime_ ) && std::isnan( ancillaryData.sequentialRangeLowestRangingComponent_ ) ) ||
        ( !std::isnan( ancillaryData.dopplerIntegrationTime_ ) && !std::isnan( ancillaryData.sequentialRangeLowestRangingComponent_ ) ) )
    {
        throw std::runtime_error(
                "Error when setting ODF metadata in tracking data: both the Doppler integration time and sequential range lowest ranging "
                "component cannot be valid values simultaneously." );
    }

    if( !std::isnan( ancillaryData.dopplerIntegrationTime_ ) )
    {
        trackingData->addAncillarySettings( "Doppler observable integration time", ancillaryData.dopplerIntegrationTime_ );
        trackingData->addAncillarySettings( "DSN Doppler reference frequency", ancillaryData.dopplerReferenceFrequency_ );
    }

    if( !std::isnan( ancillaryData.sequentialRangeLowestRangingComponent_ ) )
    {
        trackingData->addAncillarySettings( "DSN sequential range lowest ranging component",
                                            ancillaryData.sequentialRangeLowestRangingComponent_ );
    }

    trackingData->addAncillarySettings( "link ends time delays", ancillaryData.linkEndsDelays_ );
}

template< typename ObservationScalarType = double, typename TimeType = Time >
std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > > extractTrackingDataFromRawOdf(
        const std::vector< std::shared_ptr< input_output::OdfRawFileContents > >& rawOdfDataVector,
        const std::string& spacecraftName,
        const bool verboseOutput = true )
{
    auto odfFileContents = std::make_shared< PreProcessedOdfFileContents< TimeType > >( rawOdfDataVector, spacecraftName, verboseOutput );

    // Create and fill single observation sets
    std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > > trackingDataSets;

    for( auto const& [ currentObservableName, linkDataBlocks ] : odfFileContents->getProcessedDataBlocks( ) )
    {
        for( auto const& [ currentLinkEnds, currentOdfSingleLinkData ] : linkDataBlocks )
        {
            // Get vectors of times, observations, and metadata for the current observable type and link ends
            std::vector< std::vector< TimeType > > observationTimes;
            std::vector< std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > > observables;
            std::vector< OdfAncillaryData > ancillaryData;

            // Fill vectors
            separateSingleLinkOdfData( currentOdfSingleLinkData, observationTimes, observables, ancillaryData );

            // Create the single observation sets and save them
            for( unsigned int i = 0; i < observationTimes.size( ); ++i )
            {
                auto currentTrackingDataSet = std::make_shared< data::TrackingData< ObservationScalarType, TimeType > >(
                        currentObservableName, currentLinkEnds, observables.at( i ), observationTimes.at( i ), "receiver" );

                setOdfMetadataInTrackingData( currentTrackingDataSet, ancillaryData.at( i ) );
                trackingDataSets.push_back( currentTrackingDataSet );
            }
        }
    }
    return trackingDataSets;
};

template< typename ObservationScalarType = double, typename TimeType = Time >
std::vector< std::shared_ptr< data::TrackingSupplementaryData > > extractTrackingSupplementaryDataFromRawOdf(
        const std::vector< std::shared_ptr< input_output::OdfRawFileContents > >& rawOdfDataVector,
        const std::string& spacecraftName,
        const std::string& earthName,
        const bool verboseOutput = true )
{
    auto odfFileContents = std::make_shared< PreProcessedOdfFileContents< TimeType > >( rawOdfDataVector, spacecraftName, verboseOutput );

    auto const& rampTables = odfFileContents->getRampTables( );

    std::vector< std::shared_ptr< data::TrackingSupplementaryData > > trackingSupplementaryDataSets;

    for( auto const& [ stationName, rampTable ] : rampTables )
    {
        auto currentSupplementarySet = std::make_shared< data::TrackingSupplementaryData >( earthName, stationName );
        currentSupplementarySet->setFrequencySupplementaryData(
                std::vector< std::shared_ptr< data::FrequencySupplementaryData > >( { rampTable } ) );
        trackingSupplementaryDataSets.push_back( currentSupplementarySet );
    }

    return trackingSupplementaryDataSets;
};

template< typename ObservationScalarType = double, typename TimeType = Time >
std::pair< std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > >,
           std::vector< std::shared_ptr< data::TrackingSupplementaryData > > >
convertRawOdfFile( std::vector< std::shared_ptr< input_output::OdfRawFileContents > > rawOdfDataVector,
                   const std::string& spacecraftName,
                   const std::string& earthName,
                   const bool verboseOutput = true )
{
    auto const& trackingDataSets =
            extractTrackingDataFromRawOdf< ObservationScalarType, TimeType >( rawOdfDataVector, spacecraftName, verboseOutput );

    auto const& trackingSupplementaryDataSets = extractTrackingSupplementaryDataFromRawOdf< ObservationScalarType, TimeType >(
            rawOdfDataVector, spacecraftName, earthName, verboseOutput );

    return std::make_pair( trackingDataSets, trackingSupplementaryDataSets );
}

template< typename ObservationScalarType = double, typename TimeType = Time >
std::pair< std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > >,
           std::vector< std::shared_ptr< data::TrackingSupplementaryData > > >
loadOdfFile( const std::vector< std::string >& odfFileNames,
             const std::string& spacecraftName,
             const std::string& earthName,
             const bool verboseOutput = true )
{
    std::vector< std::shared_ptr< input_output::OdfRawFileContents > > rawOdfDataVector;
    for( std::string odfFileName : odfFileNames )
    {
        rawOdfDataVector.push_back( std::make_shared< input_output::OdfRawFileContents >( odfFileName ) );
    }

    return convertRawOdfFile< ObservationScalarType, TimeType >( rawOdfDataVector, spacecraftName, earthName, verboseOutput );
}

}  // namespace input_output

}  // namespace tudat

#endif  // TUDAT_PROCESSODFFILE_H

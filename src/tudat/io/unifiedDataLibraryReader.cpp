/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/unifiedDataLibraryReader.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include "tudat/astro/basic_astro/oblateSpheroidBodyShapeModel.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"

namespace tudat
{
namespace io
{

using json = nlohmann::json;

// ============================================================================
// UTASParser Implementation
// ============================================================================

template< typename ObservationScalarType, typename TimeType >
template< typename T >
T UTASParser< ObservationScalarType, TimeType >::getRequired( const json& obj, const std::string& key )
{
    if( !obj.contains( key ) )
    {
        throw std::runtime_error( "UTASParser: Required field '" + key + "' not found" );
    }

    try
    {
        return obj[ key ].get< T >( );
    }
    catch( const json::type_error& e )
    {
        throw std::runtime_error( "UTASParser: Field '" + key + "' has wrong type: " + e.what( ) );
    }
}

template< typename ObservationScalarType, typename TimeType >
template< typename T >
T UTASParser< ObservationScalarType, TimeType >::getOptional( const json& obj, const std::string& key, const T& defaultValue )
{
    if( !obj.contains( key ) )
    {
        return defaultValue;
    }

    try
    {
        return obj[ key ].get< T >( );
    }
    catch( const json::type_error& )
    {
        return defaultValue;
    }
}

template< typename ObservationScalarType, typename TimeType >
std::string UTASParser< ObservationScalarType, TimeType >::getStringOrNumber( const json& obj, const std::string& key )
{
    if( !obj.contains( key ) )
    {
        throw std::runtime_error( "UTASParser: Required field '" + key + "' not found" );
    }

    const auto& val = obj[ key ];
    if( val.is_string( ) )
    {
        return val.get< std::string >( );
    }
    else if( val.is_number_integer( ) )
    {
        return std::to_string( val.get< int64_t >( ) );
    }
    else if( val.is_number( ) )
    {
        std::ostringstream oss;
        oss << val.get< double >( );
        return oss.str( );
    }
    else
    {
        throw std::runtime_error( "UTASParser: Field '" + key + "' must be string or number" );
    }
}

template< typename ObservationScalarType, typename TimeType >
UTASParser< ObservationScalarType, TimeType >::UTASParser( const std::vector< std::string >& filePaths )
{
    this->filePaths_ = filePaths;
    parseFiles( );
}

template< typename ObservationScalarType, typename TimeType >
void UTASParser< ObservationScalarType, TimeType >::parseFiles( )
{
    for( const auto& path : this->filePaths_ )
    {
        parseFile( path );
    }

    if( observationsByStationPair_.empty( ) )
    {
        throw std::runtime_error( "UTASParser: No observations found in provided files" );
    }

    // Build concatenated vectors for legacy interface
    rebuildConcatenatedVectors( );
}

template< typename ObservationScalarType, typename TimeType >
void UTASParser< ObservationScalarType, TimeType >::rebuildConcatenatedVectors( )
{
    allEpochs_.clear( );
    allTdoa_.clear( );
    allTdoaUnc_.clear( );
    allFdoa_.clear( );
    allFdoaUnc_.clear( );

    for( const auto& entry : observationsByStationPair_ )
    {
        const auto& obs = entry.second;
        allEpochs_.insert( allEpochs_.end( ), obs.epochs.begin( ), obs.epochs.end( ) );
        allTdoa_.insert( allTdoa_.end( ), obs.tdoa.begin( ), obs.tdoa.end( ) );
        allTdoaUnc_.insert( allTdoaUnc_.end( ), obs.tdoaUnc.begin( ), obs.tdoaUnc.end( ) );
        allFdoa_.insert( allFdoa_.end( ), obs.fdoa.begin( ), obs.fdoa.end( ) );
        allFdoaUnc_.insert( allFdoaUnc_.end( ), obs.fdoaUnc.begin( ), obs.fdoaUnc.end( ) );
    }
}

template< typename ObservationScalarType, typename TimeType >
void UTASParser< ObservationScalarType, TimeType >::parseFile( const std::string& filePath )
{
    std::ifstream file( filePath );
    if( !file.is_open( ) )
    {
        throw std::runtime_error( "UTASParser: Cannot open file: " + filePath );
    }

    json j;
    try
    {
        file >> j;
    }
    catch( const json::parse_error& e )
    {
        throw std::runtime_error( "UTASParser: JSON parse error in " + filePath + ": " + e.what( ) );
    }

    // Handle different JSON structures
    json observations;
    if( j.is_array( ) )
    {
        observations = j;
    }
    else if( j.is_object( ) && j.contains( "observations" ) && j[ "observations" ].is_array( ) )
    {
        observations = j[ "observations" ];
    }
    else
    {
        throw std::runtime_error( "UTASParser: Unexpected JSON structure in " + filePath );
    }

    if( observations.empty( ) )
    {
        std::cerr << "WARNING: UTASParser: File " << filePath << " contains no observations" << std::endl;
        return;
    }

    parseObservationArray( observations, filePath );
}

template< typename ObservationScalarType, typename TimeType >
void UTASParser< ObservationScalarType, TimeType >::parseObservationArray( const json& observations, const std::string& filePath )
{
    const json& firstObs = observations[ 0 ];

    // Extract target and validate single-target constraint
    std::string targetId = getStringOrNumber( firstObs, "satNo" );
    foundTargets_.insert( targetId );
    validateSingleTarget( targetId, filePath );

    // Extract station info for this file
    std::string station1Id = getRequired< std::string >( firstObs, "origSensorId1" );
    std::string station2Id = getRequired< std::string >( firstObs, "origSensorId2" );
    StationPair stationPair = std::make_pair( station1Id, station2Id );

    // Extract station positions
    GeodeticPositionNew station1Pos;
    station1Pos.latitude = getRequired< double >( firstObs, "senlat" );
    station1Pos.longitude = getRequired< double >( firstObs, "senlon" );
    station1Pos.altitude = getRequired< double >( firstObs, "senalt" );

    GeodeticPositionNew station2Pos;
    station2Pos.latitude = getRequired< double >( firstObs, "sen2lat" );
    station2Pos.longitude = getRequired< double >( firstObs, "sen2lon" );
    station2Pos.altitude = getRequired< double >( firstObs, "sen2alt" );

    // Store station positions (will overwrite if already exists, positions should be consistent)
    stationPositions_[ station1Id ] = station1Pos;
    stationPositions_[ station2Id ] = station2Pos;

    // Initialize metadata from first file
    if( !metadataInitialized_ )
    {
        metadata_.targetId = targetId;
        metadata_.frequency = getRequired< double >( firstObs, "frequency" );
        metadata_.bandwidth = getOptional< double >( firstObs, "bandwidth", 0.0 );
        metadata_.sensor1Delay = getOptional< double >( firstObs, "sensor1Delay", 0.0 );
        metadata_.sensor2Delay = getOptional< double >( firstObs, "sensor2Delay", 0.0 );
        metadata_.dataMode = getOptional< std::string >( firstObs, "dataMode", "" );
        metadata_.origin = getOptional< std::string >( firstObs, "origin", "" );
        metadata_.source = getOptional< std::string >( firstObs, "source", "" );
        metadata_.ucts = getOptional< int >( firstObs, "ucts", 0 );

        metadataInitialized_ = true;
    }

    // Validate positions
    if( station1Pos.isZero( ) )
    {
        std::cerr << "WARNING: Station " << station1Id << " position is (0,0,0) in " << filePath << std::endl;
    }
    if( station2Pos.isZero( ) )
    {
        std::cerr << "WARNING: Station " << station2Id << " position is (0,0,0) in " << filePath << std::endl;
    }

    // Get or create observation storage for this station pair
    auto& stationObs = observationsByStationPair_[ stationPair ];

    // Reserve space for new observations
    size_t numObs = observations.size( );
    stationObs.epochs.reserve( stationObs.epochs.size( ) + numObs );
    stationObs.tdoa.reserve( stationObs.tdoa.size( ) + numObs );
    stationObs.tdoaUnc.reserve( stationObs.tdoaUnc.size( ) + numObs );
    stationObs.fdoa.reserve( stationObs.fdoa.size( ) + numObs );
    stationObs.fdoaUnc.reserve( stationObs.fdoaUnc.size( ) + numObs );

    // Parse time series data
    for( const auto& obs : observations )
    {
        // Check target consistency within file
        std::string obsTarget = getStringOrNumber( obs, "satNo" );
        foundTargets_.insert( obsTarget );
        validateSingleTarget( obsTarget, filePath );

        // Time
        std::string obTime = getRequired< std::string >( obs, "obTime" );
        TimeType epoch = convertIsoStringToEpoch( obTime );
        stationObs.epochs.push_back( epoch );

        // TDOA
        stationObs.tdoa.push_back( static_cast< ObservationScalarType >( getRequired< double >( obs, "tdoa" ) ) );
        stationObs.tdoaUnc.push_back( static_cast< ObservationScalarType >( getOptional< double >( obs, "tdoaUnc", 0.0 ) ) );

        // FDOA
        stationObs.fdoa.push_back( static_cast< ObservationScalarType >( getRequired< double >( obs, "fdoa" ) ) );
        stationObs.fdoaUnc.push_back( static_cast< ObservationScalarType >( getOptional< double >( obs, "fdoaUnc", 0.0 ) ) );
    }
}

template< typename ObservationScalarType, typename TimeType >
void UTASParser< ObservationScalarType, TimeType >::validateSingleTarget( const std::string& newTargetId, const std::string& filePath )
{
    if( foundTargets_.size( ) > 1 )
    {
        std::ostringstream oss;
        oss << "UTASParser: Multiple targets detected. BatchUTAS only supports single-target data.\n";
        oss << "Found targets: ";
        bool first = true;
        for( const auto& t : foundTargets_ )
        {
            if( !first ) oss << ", ";
            oss << "'" << t << "'";
            first = false;
        }
        oss << "\n";
        oss << "Please create separate BatchUTAS instances for each target by filtering input files.\n";
        oss << "Error occurred while parsing: " << filePath;
        throw std::runtime_error( oss.str( ) );
    }
}

template< typename ObservationScalarType, typename TimeType >
TimeType UTASParser< ObservationScalarType, TimeType >::convertIsoStringToEpoch( const std::string& isoTime )
{
    // Strip trailing 'Z' if present
    std::string timeStr = isoTime;
    if( !timeStr.empty( ) && ( timeStr.back( ) == 'Z' || timeStr.back( ) == 'z' ) )
    {
        timeStr.pop_back( );
    }

    if( timeStr.length( ) < 19 )
    {
        throw std::runtime_error( "UTASParser: Invalid time format: " + isoTime );
    }

    try
    {
        // Parse ISO string to Time
        tudat::Time timeInUTC = basic_astrodynamics::timeFromIsoString< tudat::Time >( timeStr );

        // Convert UTC to TDB
        auto timeConverter = std::make_shared< earth_orientation::TerrestrialTimeScaleConverter >( );
        Eigen::Vector3d dummyPosition( 6378.0e3, 0.0, 0.0 );
        double timeInTDB = timeConverter->getCurrentTime< double >(
            basic_astrodynamics::TimeScales::utc_scale,
            basic_astrodynamics::TimeScales::tdb_scale,
            timeInUTC,
            dummyPosition );

        return static_cast< TimeType >( timeInTDB );
    }
    catch( const std::exception& e )
    {
        throw std::runtime_error( "UTASParser: Failed to parse time '" + timeStr + "': " + e.what( ) );
    }
}

template< typename ObservationScalarType, typename TimeType >
std::vector< tudat::io::StationPair > UTASParser< ObservationScalarType, TimeType >::getStationPairs( ) const
{
    std::vector< tudat::io::StationPair > pairs;
    for( const auto& entry : observationsByStationPair_ )
    {
        pairs.push_back( entry.first );
    }
    return pairs;
}

template< typename ObservationScalarType, typename TimeType >
std::set< std::string > UTASParser< ObservationScalarType, TimeType >::getStationNames( ) const
{
    std::set< std::string > uniqueStations;
    for( const auto& entry : stationPositions_ )
    {
        uniqueStations.insert( entry.first );
    }
    return uniqueStations;
}

template< typename ObservationScalarType, typename TimeType >
Eigen::Vector3d UTASParser< ObservationScalarType, TimeType >::getStationTudatPosition( const std::string& stationName ) const
{
    auto it = stationPositions_.find( stationName );
    if( it == stationPositions_.end( ) )
    {
        throw std::runtime_error( "UTASParser: Station '" + stationName + "' not found" );
    }
    return it->second.toTudatGeodetic( );
}

template< typename ObservationScalarType, typename TimeType >
const StationPairObservations< ObservationScalarType, TimeType >&
UTASParser< ObservationScalarType, TimeType >::getObservationsForStationPair( const StationPair& stationPair ) const
{
    auto it = observationsByStationPair_.find( stationPair );
    if( it == observationsByStationPair_.end( ) )
    {
        throw std::runtime_error( "UTASParser: Station pair (" + stationPair.first + ", " + stationPair.second + ") not found" );
    }
    return it->second;
}

// ============================================================================
// BatchUTAS Implementation
// ============================================================================

template< typename ObservationScalarType, typename TimeType >
BatchUTAS< ObservationScalarType, TimeType >::BatchUTAS( const std::vector< std::string >& filePaths ):
    parser_( filePaths )
{
}

template< typename ObservationScalarType, typename TimeType >
void BatchUTAS< ObservationScalarType, TimeType >::ensureShapeModel(
    simulation_setup::SystemOfBodies& bodies,
    const std::string& stationBodyName ) const
{
    // Ensure body exists
    try
    {
        bodies.getBody( stationBodyName );
    }
    catch( const std::runtime_error& )
    {
        bodies.addBody( std::make_shared< simulation_setup::Body >( ), stationBodyName );
    }

    auto body = bodies.getBody( stationBodyName );

    // If no shape model, create oblate spheroid from SPICE
    if( body->getShapeModel( ) == nullptr )
    {
        auto shapeModel = simulation_setup::createBodyShapeModel(
            simulation_setup::fromSpiceOblateSphericalBodyShapeSettings( ), stationBodyName );
        body->setShapeModel( shapeModel );
    }
    // If shape model exists but is not oblate spheroid, throw error
    else if( std::dynamic_pointer_cast< basic_astrodynamics::OblateSpheroidBodyShapeModel >(
                 body->getShapeModel( ) ) == nullptr )
    {
        throw std::runtime_error( "BatchUTAS: Station body '" + stationBodyName +
                                  "' has incompatible shape model. Must use OblateSpheroidBodyShapeModel for ground stations." );
    }
}

template< typename ObservationScalarType, typename TimeType >
std::vector< std::string > BatchUTAS< ObservationScalarType, TimeType >::createGroundStations(
    simulation_setup::SystemOfBodies& bodies,
    const std::string& stationBodyName ) const
{
    std::vector< std::string > stationNames;
    auto body = bodies.getBody( stationBodyName );

    // Get all unique station names and their positions
    std::set< std::string > uniqueStations = parser_.getStationNames( );

    for( const auto& stationName : uniqueStations )
    {
        Eigen::Vector3d tudatPos = parser_.getStationTudatPosition( stationName );

        std::cout << "Creating ground station '" << stationName << "' on '" << stationBodyName << "'" << std::endl;
        std::cout << "    Position (alt[m], lat[rad], lon[rad]): " << tudatPos( 0 ) << ", "
                  << tudatPos( 1 ) << ", " << tudatPos( 2 ) << std::endl;

        auto settings = std::make_shared< simulation_setup::GroundStationSettings >(
            stationName, tudatPos, coordinate_conversions::geodetic_position );
        simulation_setup::createGroundStation( body, settings );
        stationNames.push_back( stationName );
    }

    return stationNames;
}

template< typename ObservationScalarType, typename TimeType >
std::vector< observation_models::LinkDefinition > BatchUTAS< ObservationScalarType, TimeType >::getLinkDefinitions(
    const std::string& stationBodyName,
    const std::string& targetNameOverride ) const
{
    const UTASMetadata& meta = parser_.getMetadata( );

    // Use override if provided, otherwise use the target ID from data
    std::string targetName = targetNameOverride.empty( ) ? meta.targetId : targetNameOverride;

    std::vector< observation_models::LinkDefinition > linkDefinitions;
    std::vector< StationPair > stationPairs = parser_.getStationPairs( );

    for( const auto& stationPair : stationPairs )
    {
        observation_models::LinkEnds linkEnds;
        linkEnds[ observation_models::receiver ] = std::make_pair( stationBodyName, stationPair.first );
        linkEnds[ observation_models::receiver2 ] = std::make_pair( stationBodyName, stationPair.second );
        linkEnds[ observation_models::transmitter ] = std::make_pair( targetName, std::string( "" ) );

        linkDefinitions.push_back( observation_models::LinkDefinition( linkEnds ) );
    }

    return linkDefinitions;
}

template< typename ObservationScalarType, typename TimeType >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > >
BatchUTAS< ObservationScalarType, TimeType >::getObservationCollection(
    const std::string& stationBodyName,
    const std::string& targetNameOverride ) const
{
    std::vector< observation_models::LinkDefinition > linkDefs = getLinkDefinitions( stationBodyName, targetNameOverride );
    std::vector< StationPair > stationPairs = parser_.getStationPairs( );

    // Create observation sets for each station pair
    std::vector< std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > > observationSetList;

    for( size_t pairIdx = 0; pairIdx < stationPairs.size( ); ++pairIdx )
    {
        const StationPair& stationPair = stationPairs[ pairIdx ];
        const observation_models::LinkDefinition& linkDef = linkDefs[ pairIdx ];

        const auto& stationObs = parser_.getObservationsForStationPair( stationPair );
        const auto& epochs = stationObs.epochs;
        const auto& tdoa = stationObs.tdoa;
        const auto& fdoa = stationObs.fdoa;

        // Build observation vectors
        std::vector< TimeType > observationTimes;
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > tdoaObservations;
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > fdoaObservations;

        observationTimes.reserve( epochs.size( ) );
        tdoaObservations.reserve( epochs.size( ) );
        fdoaObservations.reserve( epochs.size( ) );

        for( size_t i = 0; i < epochs.size( ); ++i )
        {
            observationTimes.push_back( epochs[ i ] );

            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > tdoaEntry( 1 );
            tdoaEntry( 0 ) = tdoa[ i ];
            tdoaObservations.push_back( tdoaEntry );

            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > fdoaEntry( 1 );
            fdoaEntry( 0 ) = fdoa[ i ];
            fdoaObservations.push_back( fdoaEntry );
        }

        // Create observation sets for this station pair
        auto tdoaSet = std::make_shared< observation_models::SingleObservationSet< ObservationScalarType, TimeType > >(
            observation_models::differenced_time_of_arrival,
            linkDef,
            tdoaObservations,
            observationTimes,
            observation_models::receiver );
        observationSetList.push_back( tdoaSet );

        auto fdoaSet = std::make_shared< observation_models::SingleObservationSet< ObservationScalarType, TimeType > >(
            observation_models::differenced_frequency_of_arrival,
            linkDef,
            fdoaObservations,
            observationTimes,
            observation_models::receiver );
        observationSetList.push_back( fdoaSet );
    }

    return std::make_shared< observation_models::ObservationCollection< ObservationScalarType, TimeType > >( observationSetList );
}

template< typename ObservationScalarType, typename TimeType >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > >
BatchUTAS< ObservationScalarType, TimeType >::toTudat(
    simulation_setup::SystemOfBodies& bodies,
    const std::string& stationBodyName,
    const std::string& targetNameOverride )
{
    std::string targetName = targetNameOverride.empty( ) ? getTargetId( ) : targetNameOverride;

    std::cout << "BatchUTAS: Converting to Tudat format" << std::endl;
    std::cout << "    Target ID (from data): " << getTargetId( ) << std::endl;
    if( !targetNameOverride.empty( ) )
    {
        std::cout << "    Target name (override): " << targetNameOverride << std::endl;
    }

    // Print station pairs
    std::vector< StationPair > stationPairs = parser_.getStationPairs( );
    std::cout << "    Station pairs: " << stationPairs.size( ) << std::endl;
    for( const auto& pair : stationPairs )
    {
        const auto& obs = parser_.getObservationsForStationPair( pair );
        std::cout << "        " << pair.first << " / " << pair.second
                  << " (" << obs.epochs.size( ) << " observations)" << std::endl;
    }
    std::cout << "    Total observations: " << getNumObservations( ) << std::endl;

    // Step 1: Ensure shape model
    ensureShapeModel( bodies, stationBodyName );

    // Step 2: Create ground stations
    createGroundStations( bodies, stationBodyName );

    // Step 3: Create and return observation collection
    return getObservationCollection( stationBodyName, targetNameOverride );
}

// ============================================================================
// Explicit Template Instantiations
// ============================================================================

template class UTASParser< double, double >;
template class UTASParser< double, Time >;
template class BatchUTAS< double, double >;
template class BatchUTAS< double, Time >;

}  // namespace io
}  // namespace tudat
/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_UNIFIED_DATA_LIBRARY_READER_NEW_H
#define TUDAT_UNIFIED_DATA_LIBRARY_READER_NEW_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <set>
#include <stdexcept>
#include <sstream>

#include <Eigen/Dense>
#include <nlohmann/json.hpp>

#include "tudat/basics/timeType.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/simulation/simulation.h"
#include "tudat/simulation/estimation_setup/singleObservationSet.h"
#include "tudat/simulation/estimation_setup/observationCollection.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace io
{

// ============================================================================
// GEODETIC POSITION HELPER
// ============================================================================

/**
 * @brief Geodetic position structure (altitude, latitude, longitude)
 *
 * Units depend on context - use convertToTudatGeodetic for Tudat-compatible output.
 */
struct GeodeticPositionNew
{
    double altitude;   // Length unit depends on source
    double latitude;   // Angular unit depends on source
    double longitude;  // Angular unit depends on source

    GeodeticPositionNew( ): altitude( 0.0 ), latitude( 0.0 ), longitude( 0.0 ) {}
    GeodeticPositionNew( double alt, double lat, double lon ): altitude( alt ), latitude( lat ), longitude( lon ) {}

    bool isZero( ) const
    {
        return altitude == 0.0 && latitude == 0.0 && longitude == 0.0;
    }

    Eigen::Vector3d toEigenVector( ) const
    {
        return Eigen::Vector3d( altitude, latitude, longitude );
    }

    /**
     * @brief Convert to Tudat geodetic format (altitude[m], latitude[rad], longitude[rad])
     *
     * Assumes input is in meters and degrees, converts to Tudat format.
     */
    Eigen::Vector3d toTudatGeodetic( ) const
    {
        return Eigen::Vector3d(
            altitude,
            latitude * tudat::mathematical_constants::PI / 180.0,
            longitude * tudat::mathematical_constants::PI / 180.0 );
    }
};

// ============================================================================
// UTAS METADATA STRUCT
// ============================================================================

/**
 * @brief Strongly-typed metadata for UTAS observations
 *
 * Contains target information and signal parameters common across all observations.
 * Station-specific data is stored separately as there may be multiple station pairs.
 */
struct UTASMetadata
{
    // Target identification
    std::string targetId;

    // Signal parameters (from first file - may vary slightly between files)
    double frequency = 0.0;    // Hz
    double bandwidth = 0.0;    // Hz

    // Sensor delays (from first file)
    double sensor1Delay = 0.0;  // seconds
    double sensor2Delay = 0.0;  // seconds

    // Data classification
    std::string dataMode;
    std::string origin;
    std::string source;
    int ucts = 0;

    UTASMetadata( ) = default;
};

// ============================================================================
// TYPE ALIASES
// ============================================================================

/**
 * @brief Station pair type (pair of station IDs)
 */
using StationPair = std::pair< std::string, std::string >;

/**
 * @brief Time series data for a single station pair
 */
template< typename ObservationScalarType = double, typename TimeType = double >
struct StationPairObservations
{
    std::vector< TimeType > epochs;
    std::vector< ObservationScalarType > tdoa;
    std::vector< ObservationScalarType > tdoaUnc;
    std::vector< ObservationScalarType > fdoa;
    std::vector< ObservationScalarType > fdoaUnc;

    size_t size( ) const { return epochs.size( ); }
};

// ============================================================================
// PARSER BASE CLASS
// ============================================================================

/**
 * @brief Base class for UDL format parsers
 *
 * Provides interface for parsing observation data from various formats.
 * Derived classes implement format-specific parsing logic.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class ParserBase
{
public:
    virtual ~ParserBase( ) = default;

    /**
     * @brief Parse all files and populate internal data structures
     */
    virtual void parseFiles( ) = 0;

    /**
     * @brief Get observation epochs in TDB seconds since J2000
     */
    virtual const std::vector< TimeType >& getEpochs( ) const = 0;


    /**
     * @brief Get number of observations
     */
    virtual size_t getNumObservations( ) const = 0;

    /**
     * @brief Convert geodetic position to Tudat format (radians, meters)
     *
     * @param pos Position in source units
     * @param angleInDegrees True if source angles are in degrees
     * @param altitudeInKm True if source altitude is in kilometers
     * @return Eigen::Vector3d (altitude[m], latitude[rad], longitude[rad])
     */
    static Eigen::Vector3d convertToTudatGeodetic( const GeodeticPositionNew& pos,
                                                   bool angleInDegrees = true,
                                                   bool altitudeInKm = true )
    {
        double longitude = pos.longitude;
        double latitude = pos.latitude;
        double altitude = pos.altitude;

        if( angleInDegrees )
        {
            longitude *= mathematical_constants::PI / 180.0;
            latitude *= mathematical_constants::PI / 180.0;
        }

        if( altitudeInKm )
        {
            altitude *= 1000.0;
        }

        return Eigen::Vector3d( altitude, latitude, longitude );
    }

protected:
    std::vector< std::string > filePaths_;
};

// ============================================================================
// UTAS PARSER
// ============================================================================

/**
 * @brief Parser for UTAS (Unified Tracking and Analysis System) format
 *
 * Parses JSON files in UTAS format, enforcing single-target constraint.
 * Supports multiple station pairs across files.
 * Throws std::runtime_error if multiple targets are detected.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class UTASParser : public ParserBase< ObservationScalarType, TimeType >
{
public:
    /**
     * @brief Construct parser from file paths
     *
     * @param filePaths List of JSON file paths to parse
     * @throws std::runtime_error if files contain multiple targets
     */
    explicit UTASParser( const std::vector< std::string >& filePaths );

    void parseFiles( ) override;

    // Legacy single-vector interface (concatenates all station pairs)
    const std::vector< TimeType >& getEpochs( ) const override { return allEpochs_; }
    const std::vector< ObservationScalarType >& getTdoaObservations( ) const { return allTdoa_; }
    const std::vector< ObservationScalarType >& getTdoaUncertainties( ) const { return allTdoaUnc_; }
    const std::vector< ObservationScalarType >& getFdoaObservations( ) const { return allFdoa_; }
    const std::vector< ObservationScalarType >& getFdoaUncertainties( ) const { return allFdoaUnc_; }
    size_t getNumObservations( ) const override { return allEpochs_.size( ); }

    /**
     * @brief Get parsed metadata
     */
    const UTASMetadata& getMetadata( ) const { return metadata_; }

    /**
     * @brief Get all unique station pairs
     */
    std::vector< StationPair > getStationPairs( ) const;

    /**
     * @brief Get all unique station names
     */
    std::set< std::string > getStationNames( ) const;

    /**
     * @brief Get station positions map
     */
    const std::map< std::string, GeodeticPositionNew >& getStationPositions( ) const { return stationPositions_; }

    /**
     * @brief Get observations for a specific station pair
     */
    const StationPairObservations< ObservationScalarType, TimeType >& 
    getObservationsForStationPair( const StationPair& stationPair ) const;

    /**
     * @brief Get all observations organized by station pair
     */
    const std::map< StationPair, StationPairObservations< ObservationScalarType, TimeType > >&
    getAllObservationsByStationPair( ) const { return observationsByStationPair_; }

    /**
     * @brief Get station position in Tudat format (radians, meters)
     */
    Eigen::Vector3d getStationTudatPosition( const std::string& stationName ) const;

private:
    void parseFile( const std::string& filePath );
    void parseObservationArray( const nlohmann::json& observations, const std::string& filePath );
    void validateSingleTarget( const std::string& newTargetId, const std::string& filePath );
    void rebuildConcatenatedVectors( );
    TimeType convertIsoStringToEpoch( const std::string& isoTime );

    template< typename T >
    static T getRequired( const nlohmann::json& obj, const std::string& key );

    template< typename T >
    static T getOptional( const nlohmann::json& obj, const std::string& key, const T& defaultValue );

    static std::string getStringOrNumber( const nlohmann::json& obj, const std::string& key );

    UTASMetadata metadata_;
    bool metadataInitialized_ = false;

    // Observations organized by station pair
    std::map< StationPair, StationPairObservations< ObservationScalarType, TimeType > > observationsByStationPair_;

    // Station positions (accumulated from all files)
    std::map< std::string, GeodeticPositionNew > stationPositions_;

    // Concatenated vectors for legacy interface
    std::vector< TimeType > allEpochs_;
    std::vector< ObservationScalarType > allTdoa_;
    std::vector< ObservationScalarType > allTdoaUnc_;
    std::vector< ObservationScalarType > allFdoa_;
    std::vector< ObservationScalarType > allFdoaUnc_;

    std::set< std::string > foundTargets_;  // Track all targets found for error messages
};

// ============================================================================
// BATCH UDL BASE CLASS
// ============================================================================

/**
 * @brief Base class for batch UDL observation loaders
 *
 * Provides common interface for loading observations from UDL-compliant formats.
 * Derived classes implement format-specific functionality.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class BatchUDL
{
public:
    virtual ~BatchUDL( ) = default;

    /**
     * @brief Ensure the station body has a compatible shape model
     *
     * Creates an oblate spheroid shape model from SPICE if none exists.
     * Throws if an incompatible shape model is present.
     *
     * @param bodies System of bodies
     * @param stationBodyName Name of body to check/modify
     * @throws std::runtime_error if body has incompatible shape model
     */
    virtual void ensureShapeModel( simulation_setup::SystemOfBodies& bodies,
                                   const std::string& stationBodyName ) const = 0;

    /**
     * @brief Create ground stations on the specified body
     *
     * @param bodies System of bodies (modified in place)
     * @param stationBodyName Body on which to create stations
     * @return Vector of created station names
     */
    virtual std::vector< std::string > createGroundStations(
        simulation_setup::SystemOfBodies& bodies,
        const std::string& stationBodyName ) const = 0;

    /**
     * @brief Get link definition for the observation geometry
     *
     * @param stationBodyName Name of body hosting ground stations
     * @param targetNameOverride Custom target name for link definition (empty = use ID from data)
     * @return Vector of LinkDefinitions (one per station pair)
     */
    virtual std::vector< observation_models::LinkDefinition > getLinkDefinitions(
        const std::string& stationBodyName,
        const std::string& targetNameOverride = "" ) const = 0;

    /**
     * @brief Get observation collection in Tudat format
     *
     * Creates SingleObservationSets for TDOA and FDOA observations.
     *
     * @param stationBodyName Name of body hosting ground stations
     * @param targetNameOverride Custom target name for link definition (empty = use ID from data)
     * @return Shared pointer to ObservationCollection
     */
    virtual std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > >
    getObservationCollection( const std::string& stationBodyName,
                              const std::string& targetNameOverride = "" ) const = 0;

    /**
     * @brief Full conversion to Tudat format
     *
     * Orchestrates all steps: shape model, ground stations, link definition,
     * and observation collection creation.
     *
     * @param bodies System of bodies (modified in place)
     * @param stationBodyName Body on which to create stations (default: "Earth")
     * @param targetNameOverride Custom target name for link definition (empty = use ID from data)
     * @return Shared pointer to ObservationCollection
     */
    virtual std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > >
    toTudat( simulation_setup::SystemOfBodies& bodies,
             const std::string& stationBodyName = "Earth",
             const std::string& targetNameOverride = "" ) = 0;

    /**
     * @brief Get target ID from observations
     */
    virtual std::string getTargetId( ) const = 0;

    /**
     * @brief Get station 1 ID
     */
    virtual std::string getStation1Id( ) const = 0;

    /**
     * @brief Get station 2 ID
     */
    virtual std::string getStation2Id( ) const = 0;

    /**
     * @brief Get number of observations
     */
    virtual size_t getNumObservations( ) const = 0;
};

// ============================================================================
// BATCH UTAS CLASS
// ============================================================================

/**
 * @brief Batch loader for UTAS format observations
 *
 * Main user-facing class for loading UTAS TDOA/FDOA observations.
 * Supports single-target data with multiple station pairs across files.
 * Throws if multiple targets detected.
 *
 * Example usage:
 * @code
 * BatchUTAS<double, Time> batch({"file1.json", "file2.json"});
 * auto collection = batch.toTudat(bodies, "Earth", "MySpacecraft");
 * @endcode
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class BatchUTAS : public BatchUDL< ObservationScalarType, TimeType >
{
public:
    /**
     * @brief Construct from list of JSON file paths
     *
     * @param filePaths List of UTAS JSON files
     * @throws std::runtime_error if files contain multiple targets
     */
    explicit BatchUTAS( const std::vector< std::string >& filePaths );

    // ========================================
    // BatchUDL interface implementation
    // ========================================

    void ensureShapeModel( simulation_setup::SystemOfBodies& bodies,
                           const std::string& stationBodyName ) const override;

    std::vector< std::string > createGroundStations(
        simulation_setup::SystemOfBodies& bodies,
        const std::string& stationBodyName ) const override;

    std::vector< observation_models::LinkDefinition > getLinkDefinitions(
        const std::string& stationBodyName,
        const std::string& targetNameOverride = "" ) const override;

    std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > >
    getObservationCollection( const std::string& stationBodyName,
                              const std::string& targetNameOverride = "" ) const override;

    std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > >
    toTudat( simulation_setup::SystemOfBodies& bodies,
             const std::string& stationBodyName = "Earth",
             const std::string& targetNameOverride = "" ) override;

    std::string getTargetId( ) const override { return parser_.getMetadata( ).targetId; }
    
    // These now return empty strings since there may be multiple station pairs
    std::string getStation1Id( ) const override { return ""; }
    std::string getStation2Id( ) const override { return ""; }
    
    size_t getNumObservations( ) const override { return parser_.getNumObservations( ); }

    // ========================================
    // UTAS-specific methods
    // ========================================

    /**
     * @brief Get full UTAS metadata
     */
    const UTASMetadata& getMetadata( ) const { return parser_.getMetadata( ); }

    /**
     * @brief Get all unique station pairs
     */
    std::vector< StationPair > getStationPairs( ) const { return parser_.getStationPairs( ); }

    /**
     * @brief Get all unique station names
     */
    std::set< std::string > getStationNames( ) const { return parser_.getStationNames( ); }

    /**
     * @brief Get number of station pairs
     */
    size_t getNumStationPairs( ) const { return parser_.getStationPairs( ).size( ); }

    /**
     * @brief Get observation epochs (TDB seconds since J2000) - all station pairs concatenated
     */
    const std::vector< TimeType >& getEpochs( ) const { return parser_.getEpochs( ); }

    /**
     * @brief Get TDOA observations (seconds) - all station pairs concatenated
     */
    const std::vector< ObservationScalarType >& getTdoaObservations( ) const { return parser_.getTdoaObservations( ); }

    /**
     * @brief Get TDOA uncertainties (seconds) - all station pairs concatenated
     */
    const std::vector< ObservationScalarType >& getTdoaUncertainties( ) const { return parser_.getTdoaUncertainties( ); }

    /**
     * @brief Get FDOA observations (Hz) - all station pairs concatenated
     */
    const std::vector< ObservationScalarType >& getFdoaObservations( ) const { return parser_.getFdoaObservations( ); }

    /**
     * @brief Get FDOA uncertainties (Hz) - all station pairs concatenated
     */
    const std::vector< ObservationScalarType >& getFdoaUncertainties( ) const { return parser_.getFdoaUncertainties( ); }

private:
    UTASParser< ObservationScalarType, TimeType > parser_;
};

// ============================================================================
// EXPLICIT TEMPLATE INSTANTIATION DECLARATIONS
// ============================================================================

extern template class UTASParser< double, double >;
extern template class UTASParser< double, Time >;
extern template class BatchUTAS< double, double >;
extern template class BatchUTAS< double, Time >;

}  // namespace io
}  // namespace tudat

#endif  // TUDAT_UNIFIED_DATA_LIBRARY_READER_NEW_H
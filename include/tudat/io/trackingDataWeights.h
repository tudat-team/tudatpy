/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TRACKING_DATA_WEIGHTS_H
#define TUDAT_TRACKING_DATA_WEIGHTS_H

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/io/trackingData.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace data
{

inline bool isVfcc17ValueInSet( const std::string& value, const std::set< std::string >& acceptedValues )
{
    return acceptedValues.count( value ) > 0;
}

inline std::string getVfcc17MetadataValue( const std::map< std::string, std::string >& stringMetadata, const std::string& key )
{
    const auto metadataIterator = stringMetadata.find( key );
    if( metadataIterator == stringMetadata.end( ) || metadataIterator->second.empty( ) )
    {
        return "unknown";
    }
    return metadataIterator->second;
}

inline double getVfcc17PreliminaryWeight( const double julianDate, const std::map< std::string, std::string >& stringMetadata )
{
    const std::string note2 = getVfcc17MetadataValue( stringMetadata, "note2" );
    const std::string catalog = getVfcc17MetadataValue( stringMetadata, "catalog" );
    const std::string observatory = getVfcc17MetadataValue( stringMetadata, "observatory" );

    double inverseWeight = 1000.0;

    const bool pre1890 = julianDate <= 2411368.0;
    const bool between1890And1950 = ( julianDate > 2411368.0 ) && ( julianDate <= 2433282.0 );
    const bool after1950 = julianDate > 2433282.0;

    const bool photographic = isVfcc17ValueInSet( note2, { "unknown", "P", "A", "N", "Z" } );
    if( photographic && pre1890 )
    {
        inverseWeight = 10.0;
    }
    if( photographic && between1890And1950 )
    {
        inverseWeight = 5.0;
    }
    if( photographic && after1950 )
    {
        inverseWeight = 2.5;
    }

    if( note2 == "E" || note2 == "H" )
    {
        inverseWeight = 0.2;
    }
    if( note2 == "T" )
    {
        inverseWeight = 0.5;
    }
    if( note2 == "e" )
    {
        inverseWeight = 0.75;
    }
    if( note2 == "M" )
    {
        inverseWeight = 2.0;
    }
    if( note2 == "S" || note2 == "s" )
    {
        inverseWeight = 1.5;
    }
    if( note2 == "n" )
    {
        inverseWeight = 1.0;
    }

    const bool ccd = isVfcc17ValueInSet( note2, { "C", "c", "D", "B" } );
    const bool noCatalog = catalog == "unknown";
    if( ccd && !noCatalog )
    {
        inverseWeight = 1.0;
    }
    if( ccd && noCatalog )
    {
        inverseWeight = 1.5;
    }

    if( observatory == "704" )
    {
        inverseWeight = 1.0;
    }
    if( observatory == "G96" )
    {
        inverseWeight = 0.5;
    }
    if( observatory == "F51" )
    {
        inverseWeight = 0.2;
    }
    if( observatory == "G45" )
    {
        inverseWeight = 0.6;
    }
    if( observatory == "699" )
    {
        inverseWeight = 0.8;
    }
    if( observatory == "D29" )
    {
        inverseWeight = 0.75;
    }

    if( observatory == "C51" )
    {
        inverseWeight = 1.0;
    }
    if( observatory == "E12" )
    {
        inverseWeight = 0.75;
    }
    if( observatory == "608" )
    {
        inverseWeight = 0.6;
    }
    if( observatory == "J75" )
    {
        inverseWeight = 1.0;
    }

    if( observatory == "703" )
    {
        inverseWeight = ( julianDate < 2456658.0 ) ? 1.0 : 0.8;
    }
    if( observatory == "691" )
    {
        inverseWeight = ( julianDate < 2452640.0 ) ? 0.6 : 0.5;
    }
    if( observatory == "644" )
    {
        inverseWeight = ( julianDate < 2452883.0 ) ? 0.6 : 0.4;
    }

    const bool lcoObservatory = isVfcc17ValueInSet(
            observatory, { "K92", "K93", "Q63", "Q64", "V37", "W84", "W85", "W86", "W87", "K91", "E10", "F65", "V39", "Z24", "Z31" } );
    const bool maunakeaObservatory = isVfcc17ValueInSet( observatory, { "T09", "T12", "T14" } );
    const bool ucac4Catalog = catalog == "q";
    const bool ppmxlCatalog = catalog == "t";
    const bool gaiaCatalog = isVfcc17ValueInSet( catalog, { "U", "V", "W", "X", "3", "6" } );
    const bool usnob12Catalog = isVfcc17ValueInSet( catalog, { "o", "s" } );

    if( observatory == "645" || observatory == "673" || observatory == "H01" )
    {
        inverseWeight = 0.3;
    }
    if( observatory == "689" || observatory == "950" )
    {
        inverseWeight = 0.5;
    }
    if( observatory == "J04" )
    {
        inverseWeight = 0.4;
    }
    if( observatory == "G83" && ( ucac4Catalog || ppmxlCatalog ) )
    {
        inverseWeight = 0.3;
    }
    if( observatory == "G83" && gaiaCatalog )
    {
        inverseWeight = 0.2;
    }
    if( lcoObservatory )
    {
        inverseWeight = 0.4;
    }
    if( observatory == "W84" )
    {
        inverseWeight = 0.5;
    }
    if( observatory == "Y28" && ppmxlCatalog && gaiaCatalog )
    {
        inverseWeight = 0.3;
    }
    if( observatory == "568" && usnob12Catalog )
    {
        inverseWeight = 0.5;
    }
    if( observatory == "568" && gaiaCatalog )
    {
        inverseWeight = 0.1;
    }
    if( observatory == "568" && ppmxlCatalog )
    {
        inverseWeight = 0.2;
    }
    if( maunakeaObservatory && gaiaCatalog )
    {
        inverseWeight = 0.1;
    }
    if( observatory == "309" && ( ucac4Catalog || ppmxlCatalog ) )
    {
        inverseWeight = 0.3;
    }
    if( observatory == "309" && gaiaCatalog )
    {
        inverseWeight = 0.2;
    }

    return 1.0 / std::pow( unit_conversions::convertArcSecondsToRadians( inverseWeight ), 2.0 );
}

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void setVFCC17Weights( const std::shared_ptr< TrackingData< ObservationScalarType, TimeType > > trackingData,
                       const std::vector< Eigen::Vector3d >& stationPositions,
                       const std::vector< std::map< std::string, std::string > >& stringMetadata )
{
    if( stationPositions.size( ) != trackingData->getNumberOfObservations( ) )
    {
        throw std::runtime_error( "Error when setting VFCC17 weights, number of station positions (" +
                                  std::to_string( stationPositions.size( ) ) + ") is inconsistent with number of observations (" +
                                  std::to_string( trackingData->getNumberOfObservations( ) ) + ")." );
    }
    if( stringMetadata.size( ) != trackingData->getNumberOfObservations( ) )
    {
        throw std::runtime_error( "Error when setting VFCC17 weights, number of string metadata entries (" +
                                  std::to_string( stringMetadata.size( ) ) + ") is inconsistent with number of observations (" +
                                  std::to_string( trackingData->getNumberOfObservations( ) ) + ")." );
    }

    const std::vector< TimeType > epochs = trackingData->getObservationEpochs( );
    std::vector< double > preliminaryWeights;
    preliminaryWeights.reserve( trackingData->getNumberOfObservations( ) );

    std::map< int, int > observationsPerLocalDay;
    std::vector< int > localDayIndices;
    localDayIndices.reserve( trackingData->getNumberOfObservations( ) );

    for( unsigned int i = 0; i < trackingData->getNumberOfObservations( ); ++i )
    {
        const double julianDate =
                basic_astrodynamics::convertSecondsSinceEpochToJulianDay< double >( static_cast< double >( epochs.at( i ) ) );

        double longitudeFractionOfDay = 0.0;
        if( stationPositions.at( i ).norm( ) > 0.0 )
        {
            longitudeFractionOfDay =
                    std::atan2( stationPositions.at( i ).y( ), stationPositions.at( i ).x( ) ) / ( 2.0 * mathematical_constants::PI );
        }

        const int localDayIndex = static_cast< int >( std::floor( julianDate + longitudeFractionOfDay ) );
        localDayIndices.push_back( localDayIndex );
        observationsPerLocalDay[ localDayIndex ] += 1;
        preliminaryWeights.push_back( getVfcc17PreliminaryWeight( julianDate, stringMetadata.at( i ) ) );
    }

    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > observationWeights;
    observationWeights.reserve( trackingData->getNumberOfObservations( ) );
    for( unsigned int i = 0; i < trackingData->getNumberOfObservations( ); ++i )
    {
        const double multipleObservationDeweightingFactor =
                std::max( static_cast< double >( observationsPerLocalDay.at( localDayIndices.at( i ) ) ) / 4.0, 1.0 );
        observationWeights.push_back( Eigen::Matrix< double, Eigen::Dynamic, 1 >::Constant(
                trackingData->getSingleObservationSize( ), preliminaryWeights.at( i ) / multipleObservationDeweightingFactor ) );
    }

    trackingData->setObservationWeights( observationWeights );
}

}  // namespace data

}  // namespace tudat

#endif  // TUDAT_TRACKING_DATA_WEIGHTS_H

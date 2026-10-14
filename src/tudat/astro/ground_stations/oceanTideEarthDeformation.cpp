#include <map>
#include <vector>
#include <string>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/ground_stations/oceanTideEarthDeformation.h"

namespace tudat
{

namespace basic_astrodynamics
{

//! Compatibility wrapper; neutral input records retain phases in their original degrees.
std::vector< std::map< std::string, Eigen::Matrix< double, 3, 11 > > > readBlqFile( const std::string& fileName )
{
    const auto data = input_output::readOceanLoadingData( fileName );
    std::vector< std::map< std::string, Eigen::Matrix< double, 3, 11 > > > result( 2 );
    for( const auto& station : data.stations )
    {
        result[ 0 ].emplace( station.first, station.second.amplitudesMetres );
        result[ 1 ].emplace( station.first, unit_conversions::convertDegreesToRadians( station.second.phaseLagsDegrees ) );
    }
    return result;
}

Eigen::Vector3d convertLocalOceanTideDisplacementToPlentocentricDisplacement(
        const Eigen::Vector3d& localDisplacement,
        const std::vector< Eigen::Vector3d >& earthNorthradialUnitVectors )
{
    Eigen::Vector3d planetocentricDisplacement = Eigen::Vector3d::Zero( );

    for( int i = 0; i < 3; i++ )
    {
        planetocentricDisplacement += -localDisplacement( 2 - i ) * earthNorthradialUnitVectors[ i ];
    }
    return planetocentricDisplacement;
}

//! Constructor, reads and initializes tide displacement amplitudes and phases.
OceanTideEarthDeformation::OceanTideEarthDeformation( const std::vector< std::string >& blqFiles,
                                                      const std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction ):
    doodsonArgumentFunction_( doodsonArgumentFunction )
{
    for( unsigned int i = 0; i < blqFiles.size( ); i++ )
    {
        addBlqFile( blqFiles.at( i ) );
    }

    // Initialize characteristics of 11 main tides, from Moyer(2000), Table 5-1 and Eq. 5-88
    doodsonMultipliers_ << 2.0, 0.0, 0.0, 0.0, 2.0, 2.0, -2.0, 0.0, 2.0, -1.0, 0.0, 1.0, 2.0, 2.0, 0.0, 0.0,

            1.0, 1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, 1.0, -2.0, 0.0, 1.0, -2.0, 0.0, 1.0,

            0.0, 2.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 2.0, 0.0;

    schwiderskiFactors_ << 0.0, 0.0, 0.0, 0.0, 0.5, -.50, -.50, -.50, 0.0, 0.0, 0.0;
    schwiderskiFactors_ *= mathematical_constants::PI;
}

OceanTideEarthDeformation::OceanTideEarthDeformation( const input_output::OceanLoadingData& data,
                                                      const std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction ):
    OceanTideEarthDeformation( std::vector< std::string >( ), doodsonArgumentFunction )
{
    addOceanLoadingData( data );
}

void OceanTideEarthDeformation::addBlqFile( const std::string& blqFile )
{
    addOceanLoadingData( input_output::readOceanLoadingData( blqFile ) );
}

void OceanTideEarthDeformation::addOceanLoadingData( const input_output::OceanLoadingData& data )
{
    input_output::validateOceanLoadingData( data );
    for( const auto& item : data.stations )
    {
        if( oceanLoadingData_.stations.count( item.first ) != 0 )
        {
            throw std::runtime_error( "Duplicate ocean loading station: " + item.first );
        }
    }
    oceanLoadingData_.source = oceanLoadingData_.source.empty( ) ? data.source : "multiple BLQ inputs";
    oceanLoadingData_.sourceHeaders.insert(
            oceanLoadingData_.sourceHeaders.end( ), data.sourceHeaders.begin( ), data.sourceHeaders.end( ) );
    for( const auto& item : data.stations )
    {
        oceanLoadingData_.stations.emplace( item.first, item.second );
        siteOceanTideAmplitudes_.emplace( item.first, item.second.amplitudesMetres );
        siteOceanTidePhases_.emplace( item.first, unit_conversions::convertDegreesToRadians( item.second.phaseLagsDegrees ) );
    }
}

Eigen::Vector3d OceanTideEarthDeformation::calculateDisplacementInEnuFrame( const double time, const std::string& siteIdentifier )
{
    // Check to see if requested site is available.
    if( siteOceanTideAmplitudes_.count( siteIdentifier ) == 0 )
    {
        throw std::runtime_error( "Error, requested site: " + siteIdentifier +
                                  ", which is not available for ocean tide displacement calculation (blq file missing)" );
    }

    // Calculate doodson arguments at given time.
    Eigen::Vector6d doodsonArguments = doodsonArgumentFunction_( time );

    // Initialize displacement to zero.
    Eigen::Vector3d displacement = Eigen::Vector3d::Zero( );

    Eigen::Matrix< double, 11, 1 > arguments;

    // Iterate over all tidal constituents and calculate current phase and resulting displacement.
    for( int i = 0; i < 11; i++ )
    {
        arguments[ i ] = 0.0;

        // Combine influence of doodson arguments (argument 5 and 6 are zero for all 11 constituent tides).
        for( int j = 0; j < 4; j++ )
        {
            arguments[ i ] += doodsonArguments[ j ] * doodsonMultipliers_( i, j );
        }

        // Add schwiderski phase factor for applicable tide.s
        arguments[ i ] += schwiderskiFactors_( i );

        // Calculate site displacement due to tide in three dimensions.
        for( int j = 0; j < 3; j++ )
        {
            displacement( j ) += siteOceanTideAmplitudes_[ siteIdentifier ]( j, i ) *
                    std::cos( arguments[ i ] - siteOceanTidePhases_[ siteIdentifier ]( j, i ) );
        }
    }

    return ( Eigen::Vector3d( ) << -displacement( 1 ), -displacement( 2 ), displacement( 0 ) ).finished( );
}

//! Function to calculate site displacement at given time.
Eigen::Vector3d OceanTideEarthDeformation::calculateDisplacement(
        const double time,
        const std::shared_ptr< ground_stations::GroundStationState > stationState )
{
    return stationState->getRotationFromBodyFixedToTopocentricFrame( time ).inverse( ) *
            calculateDisplacementInEnuFrame( time, stationState->getSiteId( ) );
}

}  // namespace basic_astrodynamics

}  // namespace tudat

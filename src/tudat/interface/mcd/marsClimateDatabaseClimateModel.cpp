#include "tudat/interface/mcd/marsClimateDatabaseClimateModel.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/basic_astro/timeConversions.h"

namespace tudat
{

namespace mcd_interface
{

MarsClimateDatabaseClimateModel::MarsClimateDatabaseClimateModel( const std::string& mcdDataPath,
                                                                  const int dustScenario,
                                                                  const int perturbationKey,
                                                                  const double perturbationSeed,
                                                                  const double gravityWaveLength,
                                                                  const int highResolutionMode,
                                                                  const int maximumCacheSize ):
    ClimateModel( ), mcdDataPath_( mcdDataPath ), dustScenario_( dustScenario ), perturbationKey_( perturbationKey ),
    perturbationSeed_( perturbationSeed ), gravityWaveLength_( gravityWaveLength ), highResolutionMode_( highResolutionMode ),
    maximumCacheSize_( maximumCacheSize ), density_( 0.0 ), pressure_( 0.0 ), temperature_( 0.0 ), zonalWind_( 0.0 ), meridionalWind_( 0.0 )
{
    zkey_ = 3;

    // Set default MCD data path if not provided
    if( mcdDataPath_.empty( ) )
    {
        // Use the compile-time defined path to the MCD data directory
        // This is set by CMake and points to the source tree location
#ifdef MCD_DATA_PATH
        mcdDataPath_ = MCD_DATA_PATH;
#else
        // Fallback: this should not normally happen if CMake is configured correctly
        throw std::runtime_error(
                "Error in MCD climate model: MCD_DATA_PATH not defined at compile time. "
                "Please ensure MCD is properly configured in CMake." );
#endif
    }

    // Ensure path ends with '/'
    if( !mcdDataPath_.empty( ) && mcdDataPath_.back( ) != '/' )
    {
        mcdDataPath_ += '/';
    }

    if( ( dustScenario_ < 1 || dustScenario_ > 8 ) && ( dustScenario_ < 24 || dustScenario_ > 35 ) )
    {
        throw std::runtime_error( "McdClimateModel: Invalid dustScenario " + std::to_string( dustScenario_ ) + ". Must be 1-8 or 24-35." );
    }

    if( perturbationKey_ < 0 || perturbationKey_ > 5 )
    {
        throw std::runtime_error( "McdClimateModel: Invalid perturbationKey " + std::to_string( perturbationKey_ ) + ". Must be 0-5." );
    }

    // Validate perturbationSeed for perturbationKey=5
    if( perturbationKey_ == 5 )
    {
        if( perturbationSeed_ < -4.0 || perturbationSeed_ > 4.0 )
        {
            throw std::runtime_error(
                    "McdClimateModel: For perturbationKey=5, perturbationSeed must be in [-4, 4]. "
                    "Got: " +
                    std::to_string( perturbationSeed_ ) );
        }
    }

    if( highResolutionMode_ != 0 && highResolutionMode_ != 1 )
    {
        throw std::runtime_error( "McdClimateModel: Invalid highResolutionMode " + std::to_string( highResolutionMode_ ) +
                                  ". Must be 0 or 1." );
    }

    if( maximumCacheSize_ < 1 )
    {
        throw std::runtime_error( "McdClimateModel: Invalid maximumCacheSize " + std::to_string( maximumCacheSize_ ) +
                                  ". Must be at least 1." );
    }

    // Validate gravity wave length if perturbations are used
    if( ( perturbationKey_ == 3 || perturbationKey_ == 4 ) && gravityWaveLength_ < 0.0 )
    {
        throw std::runtime_error(
                "McdClimateModel: gravityWaveLength must be >= 0.0 when using gravity wave perturbations. "
                "Got: " +
                std::to_string( gravityWaveLength_ ) );
    }
}

void MarsClimateDatabaseClimateModel::setZkey( int zkey )
{
    if( zkey < 1 || zkey > 4 )
    {
        throw std::runtime_error( "McdClimateModel: Invalid zkey " + std::to_string( zkey ) + ". Must be 1-4." );
    }
    if( zkey_ != zkey )
    {
        zkey_ = zkey;
        mcdCache_.clear( );
        mcdCacheInsertionOrder_.clear( );
    }
}

void MarsClimateDatabaseClimateModel::updateCache( const double verticalCoordinate,
                                                   const double longitude,
                                                   const double latitude,
                                                   const double time )
{
    const McdCacheKey cacheKey = getCacheKey( verticalCoordinate, longitude, latitude, time );
    if( mcdCache_.count( cacheKey ) > 0 )
    {
        return;
    }

    int dateKey = 0;
    double xdate = basic_astrodynamics::convertSecondsSinceEpochToJulianDay< double >( time );
    float localTime = 0.0f;
    float seedin_f = static_cast< float >( perturbationSeed_ );
    float gwlength_f = static_cast< float >( gravityWaveLength_ );
    float verticalCoordinate_f = static_cast< float >( verticalCoordinate );
    float longitudeDeg = static_cast< float >( unit_conversions::convertRadiansToDegrees( longitude ) );
    float latitudeDeg = static_cast< float >( unit_conversions::convertRadiansToDegrees( latitude ) );
    float seedout;
    int ier;

    TUDAT_MCD_FORTRAN_CALL_MCD( &zkey_,
                                &verticalCoordinate_f,
                                &longitudeDeg,
                                &latitudeDeg,
                                &highResolutionMode_,
                                &dateKey,
                                &xdate,
                                &localTime,
                                mcdDataPath_.c_str( ),
                                &dustScenario_,
                                &perturbationKey_,
                                &seedin_f,
                                &gwlength_f,
                                extraVariableKeys_,
                                &pressure_,
                                &density_,
                                &temperature_,
                                &zonalWind_,
                                &meridionalWind_,
                                meanVariables_,
                                extraVariables_,
                                &seedout,
                                &ier,
                                mcdDataPath_.length( ) );

    if( ier != 0 )
    {
        throw std::runtime_error( "McdAtmosphereModel: MCD routine returned error code " + std::to_string( ier ) );
    }

    std::vector< double > meanVariables( 5 );
    std::vector< double > extraVariables( 100 );
    for( int i = 0; i < 5; i++ )
    {
        meanVariables[ i ] = static_cast< double >( meanVariables_[ i ] );
    }
    for( int i = 0; i < 100; i++ )
    {
        extraVariables[ i ] = static_cast< double >( extraVariables_[ i ] );
    }

    mcdCache_[ cacheKey ] = std::make_shared< McdCache >( density_, pressure_, temperature_, zonalWind_, meridionalWind_ );
    mcdCache_[ cacheKey ]->meanVariables_ = meanVariables;
    mcdCache_[ cacheKey ]->extraVariables_ = extraVariables;
    mcdCacheInsertionOrder_.push_back( cacheKey );
    while( static_cast< int >( mcdCache_.size( ) ) > maximumCacheSize_ )
    {
        mcdCache_.erase( mcdCacheInsertionOrder_.front( ) );
        mcdCacheInsertionOrder_.pop_front( );
    }
}

void MarsClimateDatabaseClimateModel::addExtraVariableKeys( std::vector< mcd_interface::ExtVar > requiredExtraVariables )
{
    for( auto it : requiredExtraVariables )
    {
        extraVariableKeys_[ static_cast< int >( it ) ] = 1;
    }
    mcdCache_.clear( );
    mcdCacheInsertionOrder_.clear( );
}

}  // namespace mcd_interface

}  // namespace tudat

#include "tudat/astro/ground_stations/stationLoadingDeformation.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include <algorithm>
#include <limits>

namespace tudat
{
namespace basic_astrodynamics
{

StationLoadingDeformation::StationLoadingDeformation( const input_output::StationLoadingData& data,
                                                      const std::map< std::string, std::string >& stationBindings,
                                                      const double maximumGapSeconds,
                                                      std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter ):
    BodyDeformationModel( converter ), data_( data ), stationBindings_( stationBindings ), maximumGapSeconds_( maximumGapSeconds )
{
    for( const auto& binding : stationBindings_ )
    {
        if( !data_.sites.count( binding.second ) )
        {
            throw std::invalid_argument( "Station loading binding names an unknown site " + binding.second );
        }
    }
    for( const auto& tai : data_.epochs )
    {
        tdbEpochs_.push_back( timeConverter_->getCurrentTime< Time >( tai_scale, tdb_scale, tai ) );
    }
    for( const auto& entry : data_.sites )
    {
        const Eigen::Vector3d& position = entry.second.cartesianPosition;
        Eigen::Matrix3d basis;
        basis.col( 0 ) = position.normalized( );
        basis.col( 1 ) = Eigen::Vector3d( -position.y( ), position.x( ), 0 ).normalized( );
        basis.col( 2 ) = basis.col( 0 ).cross( basis.col( 1 ) );
        siteBases_.emplace( entry.first, basis );
    }
}

Eigen::Vector6d StationLoadingDeformation::evaluate( const Time& requestedTdb,
                                                     const std::shared_ptr< ground_stations::GroundStationState >& station,
                                                     const bool withRate ) const
{
    // Select the site bound to the station, or the nearest site within whose radius the station lies.
    const Eigen::Vector3d position = station->getNominalCartesianPosition( );
    const auto binding = stationBindings_.find( station->getSiteId( ) );
    if( !stationBindings_.empty( ) && binding == stationBindings_.end( ) )
    {
        throw std::out_of_range( "No station loading binding for station " + station->getSiteId( ) );
    }
    std::string identifier;
    double distance = std::numeric_limits< double >::infinity( );
    for( const auto& entry : data_.sites )
    {
        const double current = ( position - entry.second.cartesianPosition ).norm( );
        if( current <= entry.second.validityRadiusMetres && current < distance &&
            ( binding == stationBindings_.end( ) || binding->second == entry.first ) )
        {
            identifier = entry.first;
            distance = current;
        }
    }
    if( identifier.empty( ) )
    {
        throw std::out_of_range( "No station loading site covers station " + station->getSiteId( ) );
    }
    const auto& displacements = data_.sites.at( identifier ).displacements;
    if( displacements.empty( ) )
    {
        throw std::out_of_range( "Station loading site " + identifier + " has no displacements" );
    }

    // An epoch less than one microsecond outside the first or last sample is evaluated at that sample.
    const std::size_t first = displacements.begin( )->first, last = displacements.rbegin( )->first;
    if( requestedTdb < tdbEpochs_[ first ] - 1.0E-6 || requestedTdb > tdbEpochs_[ last ] + 1.0E-6 )
    {
        throw std::out_of_range( "Station loading epoch outside the series of site " + identifier );
    }
    const Time tdb = std::min( std::max( requestedTdb, tdbEpochs_[ first ] ), tdbEpochs_[ last ] );
    const auto node = std::lower_bound( tdbEpochs_.begin( ) + first, tdbEpochs_.begin( ) + last + 1, tdb );
    const std::size_t index = std::distance( tdbEpochs_.begin( ), node );
    const bool atSample = *node == tdb;

    // Interpolate in TAI, using the TAI epoch of the series at a sample.
    Eigen::Vector6d result = Eigen::Vector6d::Zero( );
    if( atSample )
    {
        result.head< 3 >( ) = displacements.at( index );
    }
    if( !atSample || withRate )
    {
        if( first == last )
        {
            throw std::out_of_range( "Station loading rate unavailable for a single sample" );
        }
        const std::size_t lower = atSample ? ( index == last ? last - 1 : index ) : index - 1;
        const double duration = static_cast< double >( data_.epochs[ lower + 1 ] - data_.epochs[ lower ] );
        if( duration > maximumGapSeconds_ )
        {
            throw std::out_of_range( "Station loading interval exceeds the maximum gap" );
        }
        const Eigen::Vector3d slope = ( displacements.at( lower + 1 ) - displacements.at( lower ) ) / duration;
        if( !atSample )
        {
            const Time tai = timeConverter_->getTaiTime( tdb_scale, tdb );
            result.head< 3 >( ) = displacements.at( lower ) + slope * static_cast< double >( tai - data_.epochs[ lower ] );
        }
        if( withRate )
        {
            result.tail< 3 >( ) = slope * timeConverter_->getGeocentricTtRateWithRespectToTdb( tdb );
        }
    }
    const auto& basis = siteBases_.at( identifier );
    result.head< 3 >( ) = basis * result.head< 3 >( ).eval( );
    result.tail< 3 >( ) = basis * result.tail< 3 >( ).eval( );
    return result;
}

Eigen::Vector3d StationLoadingDeformation::evaluateDisplacement( const Time& tdb,
                                                                 const std::shared_ptr< ground_stations::GroundStationState > station )
{
    return evaluate( tdb, station, false ).head< 3 >( );
}

Eigen::Vector6d StationLoadingDeformation::evaluateDisplacementState( const Time& tdb,
                                                                      const std::shared_ptr< ground_stations::GroundStationState > station )
{
    return evaluate( tdb, station, true );
}

}  // namespace basic_astrodynamics
}  // namespace tudat

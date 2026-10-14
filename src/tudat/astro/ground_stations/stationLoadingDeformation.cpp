#include "tudat/astro/ground_stations/stationLoadingDeformation.h"
#include "tudat/astro/ground_stations/atmosphericTideDeformation.h"
#include "tudat/astro/ground_stations/oceanTideEarthDeformation.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include <algorithm>
#include <cmath>
#include <set>
#include <limits>

namespace tudat
{
namespace basic_astrodynamics
{
namespace
{
using namespace input_output;
Eigen::Vector3d vector( const std::array< double, 3 >& values )
{
    return { values[ 0 ], values[ 1 ], values[ 2 ] };
}
bool hasComponent( const StationLoadingMetadata& metadata, const StationLoadingComponent component )
{
    return std::find( metadata.components.begin( ), metadata.components.end( ), component ) != metadata.components.end( );
}
bool sameScope( const StationLoadingMetadata& first, const StationLoadingMetadata& second )
{
    if( first.components.empty( ) ) return false;
    for( const auto signal : { first.diurnal, first.semidiurnal, first.annual, first.semiannual, first.mean, first.trend } )
        if( signal == StationLoadingSignalPresence::unknown ) return false;
    return std::set< StationLoadingComponent >( first.components.begin( ), first.components.end( ) ) ==
            std::set< StationLoadingComponent >( second.components.begin( ), second.components.end( ) ) &&
            first.diurnal == second.diurnal && first.semidiurnal == second.semidiurnal && first.annual == second.annual &&
            first.semiannual == second.semiannual && first.mean == second.mean && first.trend == second.trend &&
            std::set< std::string >( first.removedHarmonics.begin( ), first.removedHarmonics.end( ) ) ==
            std::set< std::string >( second.removedHarmonics.begin( ), second.removedHarmonics.end( ) ) &&
            first.oceanResponse == second.oceanResponse;
}
}  // namespace

StationLoadingDeformation::StationLoadingDeformation( const input_output::StationLoadingData& data,
                                                      const StationLoadingEvaluationSettings& settings,
                                                      std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter ):
    BodyDeformationModel( converter ), data_( data ), settings_( settings )
{
    using namespace input_output;
    validateStationLoadingData( data_ );
    if( !std::isfinite( settings_.maximumGapSeconds ) || settings_.maximumGapSeconds <= 0 )
        throw std::invalid_argument( "Loading maximum interpolation gap must be positive and finite" );
    if( settings_.spatialSelection != StationLoadingSpatialSelection::explicitBindings &&
        settings_.spatialSelection != StationLoadingSpatialSelection::withinProviderRadius )
        throw std::invalid_argument( "Unknown loading spatial-selection policy" );
    if( settings_.spatialSelection == StationLoadingSpatialSelection::explicitBindings && settings_.stationBindings.empty( ) )
        throw std::invalid_argument( "Loading requires explicit station/provider bindings or requested spatial selection" );
    if( settings_.spatialSelection == StationLoadingSpatialSelection::withinProviderRadius && !settings_.stationBindings.empty( ) )
        throw std::invalid_argument( "Spatial loading selection cannot silently combine with explicit bindings" );
    for( const auto& binding : settings_.stationBindings )
        if( !data_.sites.count( binding.second ) ) throw std::invalid_argument( "Loading binding names an undeclared provider label" );
    if( data_.basis == StationLoadingBasis::unspecified || data_.metadata.referenceFrame.empty( ) ||
        data_.metadata.referenceOrigin == StationLoadingOrigin::unspecified || data_.metadata.originDefinition.empty( ) ||
        settings_.outputOriginDefinition.empty( ) )
        throw std::invalid_argument(
                "Loading basis, Cartesian reference frame, source origin and requested output origin must be explicit" );
    if( !settings_.allowUnknownInventory )
    {
        if( data_.metadata.components.empty( ) || data_.metadata.product.empty( ) || data_.metadata.release.empty( ) )
            throw std::invalid_argument(
                    "Loading source product/release/component inventory is unknown; resolve it or explicitly accept unknown inventory" );
        for( const auto signal : { data_.metadata.diurnal,
                                   data_.metadata.semidiurnal,
                                   data_.metadata.annual,
                                   data_.metadata.semiannual,
                                   data_.metadata.mean,
                                   data_.metadata.trend } )
            if( signal == StationLoadingSignalPresence::unknown )
                throw std::invalid_argument( "Unknown loading signal inventory requires an explicit decision" );
    }
    grid_ = prepareGrid( data_ );
    for( const auto& entry : data_.sites )
    {
        Eigen::Matrix3d basis = Eigen::Matrix3d::Identity( );
        if( data_.basis == StationLoadingBasis::geocentricRadialEastNorth )
        {
            const auto position = vector( entry.second.cartesianPosition );
            const double norm = position.norm( ), horizontal = std::hypot( position.x( ), position.y( ) );
            if( norm == 0 || horizontal == 0 )
                throw std::invalid_argument( "Provider radial/east/north axes are undefined at geocentre/pole" );
            basis.col( 0 ) = position / norm;
            basis.col( 1 ) = Eigen::Vector3d( -position.y( ), position.x( ), 0 ) / horizontal;
            basis.col( 2 ) = basis.col( 0 ).cross( basis.col( 1 ) );
        }
        sourceBases_.emplace( entry.first, basis );
    }
    if( settings_.originTranslation )
    {
        // Freeze preloaded translation data too; later caller mutation cannot change this model.
        settings_.originTranslation = std::make_shared< const StationLoadingOriginTranslation >( *settings_.originTranslation );
        const auto& translation = *settings_.originTranslation;
        validateStationLoadingData( translation.series );
        if( translation.identifier.empty( ) || translation.sourceOriginDefinition != data_.metadata.originDefinition ||
            translation.targetOriginDefinition != settings_.outputOriginDefinition ||
            translation.targetOrigin == StationLoadingOrigin::unspecified ||
            translation.series.metadata.originDefinition != translation.targetOriginDefinition ||
            translation.series.metadata.referenceOrigin != translation.targetOrigin ||
            translation.series.metadata.referenceFrame != data_.metadata.referenceFrame ||
            translation.series.basis != StationLoadingBasis::cartesian || translation.series.sites.size( ) != 1 ||
            !sameScope( data_.metadata, translation.series.metadata ) )
            throw std::invalid_argument(
                    "Loading origin translation must match both definitions, Cartesian axes and the complete component/signal scope" );
        translationGrid_ = prepareGrid( translation.series );
    }
    else if( settings_.outputOriginDefinition != data_.metadata.originDefinition )
        throw std::invalid_argument( "Loading origin change requires a matching explicit Cartesian translation series" );
}

StationLoadingDeformation::PreparedGrid StationLoadingDeformation::prepareGrid( const input_output::StationLoadingData& data ) const
{
    PreparedGrid result;
    for( const auto& source : data.epochs )
    {
        const long long hours = ( static_cast< long long >( source.modifiedJulianDay ) - 51544LL ) * 24LL - 12LL;
        if( hours < std::numeric_limits< int >::min( ) || hours > std::numeric_limits< int >::max( ) )
            throw std::invalid_argument( "Loading source epoch is outside the precise Time range" );
        const Time time( static_cast< int >( hours ), source.secondsOfDay );
        Time tai;
        switch( data.timeScale )
        {
            case input_output::StationLoadingTimeScale::tai:
                tai = time;
                break;
            case input_output::StationLoadingTimeScale::tt:
                tai = convertTTtoTAI( time );
                break;
            case input_output::StationLoadingTimeScale::utc:
                tai = sofa_interface::convertUTCtoTAI( time );
                break;
            case input_output::StationLoadingTimeScale::tdb:
                tai = timeConverter_->getTaiTime( tdb_scale, time );
                break;
            default:
                throw std::invalid_argument( "Unsupported source loading clock" );
        }
        if( !result.tai.empty( ) && tai <= result.tai.back( ) )
            throw std::invalid_argument( "Source loading epochs are not increasing on continuous TAI" );
        result.tai.push_back( tai );
        result.tdb.push_back( timeConverter_->getCurrentTime< Time >( tai_scale, tdb_scale, tai ) );
    }
    return result;
}

const input_output::StationLoadingSite& StationLoadingDeformation::selectSite( const Eigen::Vector3d& position,
                                                                               const std::string& stationId ) const
{
    if( !position.allFinite( ) ) throw std::invalid_argument( "Loading target coordinates must be finite" );
    const input_output::StationLoadingSite* selected = nullptr;
    const auto inRadius = [ & ]( const auto& site ) {
        return ( position - vector( site.cartesianPosition ) ).norm( ) <= site.validityRadiusMetres;
    };
    if( settings_.spatialSelection == StationLoadingSpatialSelection::explicitBindings )
    {
        const auto binding = settings_.stationBindings.find( stationId );
        if( binding == settings_.stationBindings.end( ) )
            throw std::out_of_range( "No explicit loading binding for this environment station ID" );
        selected = &data_.sites.at( binding->second );
        if( !inRadius( *selected ) ) throw std::out_of_range( "Loading target is outside the bound provider validity radius" );
    }
    else
        for( const auto& entry : data_.sites )
            if( inRadius( entry.second ) )
            {
                if( selected ) throw std::invalid_argument( "Loading spatial match is ambiguous within provider validity radii" );
                selected = &entry.second;
            }
    if( !selected ) throw std::out_of_range( "No loading provider point covers this target" );
    return *selected;
}

Eigen::Vector3d StationLoadingDeformation::interpolate( const input_output::StationLoadingSite& site,
                                                        const PreparedGrid& grid,
                                                        const Time& tdb,
                                                        const bool rate ) const
{
    if( site.samples.empty( ) ) throw std::out_of_range( "Loading site has no source samples" );
    const std::size_t first = site.samples.begin( )->first, last = site.samples.rbegin( )->first;
    if( tdb < grid.tdb[ first ] || tdb > grid.tdb[ last ] ) throw std::out_of_range( "Loading epoch outside exact site coverage" );
    const auto sourceNode = std::lower_bound( grid.tdb.begin( ) + first, grid.tdb.begin( ) + last + 1, tdb );
    const bool exactNode = *sourceNode == tdb;
    const std::size_t sourceIndex = std::distance( grid.tdb.begin( ), sourceNode );
    const Time tai = exactNode ? grid.tai[ sourceIndex ] : timeConverter_->getTaiTime( tdb_scale, tdb );
    if( exactNode && !rate ) return vector( site.samples.at( sourceIndex ).displacement );
    if( first == last ) throw std::out_of_range( "Loading position exists at this singleton epoch, but no interpolant rate is available" );
    const std::size_t lower = exactNode ? ( sourceIndex == last ? last - 1 : sourceIndex ) : sourceIndex - 1;
    const double duration = static_cast< double >( grid.tai[ lower + 1 ] - grid.tai[ lower ] );
    if( duration > settings_.maximumGapSeconds ) throw std::out_of_range( "Loading interpolant support exceeds the declared maximum gap" );
    const Eigen::Vector3d left = vector( site.samples.at( lower ).displacement );
    const Eigen::Vector3d right = vector( site.samples.at( lower + 1 ).displacement );
    if( rate ) return ( right - left ) / duration * timeConverter_->getGeocentricTtRateWithRespectToTdb( tdb );
    return left + ( right - left ) * ( static_cast< double >( tai - grid.tai[ lower ] ) / duration );
}

Eigen::Vector3d StationLoadingDeformation::evaluate( const Time& tdb,
                                                     const Eigen::Vector3d& position,
                                                     const std::string& stationId,
                                                     const bool rate ) const
{
    const auto& site = selectSite( position, stationId );
    Eigen::Vector3d result = sourceBases_.at( site.identifier ) * interpolate( site, grid_, tdb, rate );
    if( settings_.originTranslation )
        result += interpolate( settings_.originTranslation->series.sites.begin( )->second, translationGrid_, tdb, rate );
    return result;
}

Eigen::Vector3d StationLoadingDeformation::evaluateDisplacement( const Time& tdb,
                                                                 const std::shared_ptr< ground_stations::GroundStationState > station )
{
    return evaluate( tdb, station->getNominalCartesianPosition( ), station->getSiteId( ), false );
}
Eigen::Vector6d StationLoadingDeformation::evaluateDisplacementState( const Time& tdb,
                                                                      const std::shared_ptr< ground_stations::GroundStationState > station )
{
    Eigen::Vector6d state;
    state.head< 3 >( ) = evaluateDisplacement( tdb, station );
    state.tail< 3 >( ) = evaluate( tdb, station->getNominalCartesianPosition( ), station->getSiteId( ), true );
    return state;
}
Time StationLoadingDeformation::normalizeScalarTdb( const double tdb ) const
{
    std::optional< Time > exact;
    for( const auto* grid : { &grid_, &translationGrid_ } )
    {
        const auto found = std::lower_bound( grid->tdb.begin( ), grid->tdb.end( ), tdb, []( const Time& source, double input ) {
            return static_cast< double >( source ) < input;
        } );
        if( found == grid->tdb.end( ) || static_cast< double >( *found ) != tdb ) continue;
        if( ( found + 1 != grid->tdb.end( ) && static_cast< double >( *( found + 1 ) ) == tdb ) || ( exact && *exact != *found ) )
            throw std::invalid_argument( "Scalar TDB input cannot distinguish loading source epochs; use precise Time input" );
        exact = *found;
    }
    return exact ? *exact : BodyDeformationModel::normalizeScalarTdb( tdb );
}

void validateStationLoadingComposition( const std::vector< std::shared_ptr< BodyDeformationModel > >& models,
                                        const std::vector< bool >& selected,
                                        const std::shared_ptr< ground_stations::GroundStationState >& station )
{
    using namespace input_output;
    if( models.size( ) != selected.size( ) ) throw std::invalid_argument( "Loading composition selection count mismatch" );
    std::set< std::string > translations;
    for( std::size_t i = 0; i < models.size( ); ++i )
    {
        const auto loading = std::dynamic_pointer_cast< StationLoadingDeformation >( models[ i ] );
        if( !selected[ i ] || !loading ) continue;
        const auto& metadata = loading->getData( ).metadata;
        const auto& settings = loading->getEvaluationSettings( );
        if( settings.originTranslation && !translations.insert( settings.originTranslation->identifier ).second )
            throw std::invalid_argument( "The same global loading origin translation would be applied more than once" );
        for( std::size_t j = i + 1; j < models.size( ); ++j )
            if( selected[ j ] )
                if( const auto other = std::dynamic_pointer_cast< StationLoadingDeformation >( models[ j ] ) )
                    if( settings.outputOriginDefinition != other->getEvaluationSettings( ).outputOriginDefinition ||
                        metadata.referenceFrame != other->getData( ).metadata.referenceFrame )
                        throw std::invalid_argument(
                                "Selected loading contributions do not share explicit Cartesian axes and output origin definitions" );
        if( settings.allowDeclaredOverlap ) continue;
        if( station && station->hasFittedSeasonalTerms( ) &&
            ( metadata.annual != StationLoadingSignalPresence::absent || metadata.semiannual != StationLoadingSignalPresence::absent ) )
            throw std::invalid_argument(
                    "Physical seasonal loading and fitted station seasonal terms require an explicit overlap decision" );
        for( std::size_t j = 0; j < models.size( ); ++j )
        {
            if( i == j || !selected[ j ] ) continue;
            if( const auto other = std::dynamic_pointer_cast< StationLoadingDeformation >( models[ j ] ) )
            {
                for( const auto component : metadata.components )
                    if( hasComponent( other->getData( ).metadata, component ) && !other->getEvaluationSettings( ).allowDeclaredOverlap )
                        throw std::invalid_argument( "Selected loading series declare the same physical component twice" );
            }
            const bool atmospheric = hasComponent( metadata, StationLoadingComponent::atmosphere ) &&
                    std::dynamic_pointer_cast< AtmosphericTideDeformation >( models[ j ] ) != nullptr;
            const bool ocean = hasComponent( metadata, StationLoadingComponent::ocean ) &&
                    std::dynamic_pointer_cast< OceanTideEarthDeformation >( models[ j ] ) != nullptr;
            const bool retainedTides = metadata.diurnal == StationLoadingSignalPresence::present ||
                    metadata.semidiurnal == StationLoadingSignalPresence::present;
            const bool unknownTides = metadata.diurnal == StationLoadingSignalPresence::unknown ||
                    metadata.semidiurnal == StationLoadingSignalPresence::unknown;
            if( ( atmospheric || ocean ) && ( retainedTides || ( unknownTides && !settings.allowUnknownInventory ) ) )
                throw std::invalid_argument( "Loading series retains tidal bands already supplied by another selected loading model" );
        }
    }
}
}  // namespace basic_astrodynamics
}  // namespace tudat

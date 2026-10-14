/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/astro/ground_stations/stationReferencePoint.h"
#include "tudat/astro/ground_stations/transformedStationHistoryTrajectory.h"
#include "tudat/astro/basic_astro/geodeticCoordinateConversions.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace tudat
{
namespace ground_stations
{
namespace
{
using namespace input_output;
void require( bool condition, const std::string& message )
{
    if( !condition ) throw std::runtime_error( "Station reference point: " + message );
}
bool sameSite( const SinexSiteKey& a, const SinexSiteKey& b )
{
    return !( a < b ) && !( b < a );
}
bool sameIdentity( const StationReferencePointIdentity& a, const StationReferencePointIdentity& b )
{
    return sameSite( a.site, b.site ) && a.domesId == b.domesId && a.kind == b.kind && a.name == b.name;
}
bool contains( const StationHistoryInterval& interval, const Time& tai )
{
    return ( !interval.startTai || tai >= *interval.startTai ) && ( !interval.endTai || tai < *interval.endTai );
}
bool overlap( const StationHistoryInterval& a, const StationHistoryInterval& b )
{
    return !( a.endTai && b.startTai && *a.endTai <= *b.startTai ) && !( b.endTai && a.startTai && *b.endTai <= *a.startTai );
}
StationHistoryInterval intersection( const StationHistoryInterval& a, const StationHistoryInterval& b )
{
    return { !a.startTai ? b.startTai : ( !b.startTai ? a.startTai : std::optional< Time >( std::max( *a.startTai, *b.startTai ) ) ),
             !a.endTai ? b.endTai : ( !b.endTai ? a.endTai : std::optional< Time >( std::min( *a.endTai, *b.endTai ) ) ) };
}
template< typename Record >
void sortAndCheck( std::vector< Record >& records )
{
    require( !records.empty( ), "empty offset/equipment history" );
    std::sort( records.begin( ), records.end( ), []( const auto& a, const auto& b ) {
        return !a.validity.startTai ? b.validity.startTai.has_value( )
                                    : ( b.validity.startTai && *a.validity.startTai < *b.validity.startTai );
    } );
    for( std::size_t i = 0; i < records.size( ); ++i )
    {
        const auto& interval = records[ i ].validity;
        for( const auto& boundary : { interval.startTai, interval.endTai } )
            require( !boundary || std::isfinite( static_cast< long double >( *boundary ) ), "nonfinite validity endpoint" );
        require( !interval.startTai || !interval.endTai || *interval.startTai < *interval.endTai, "empty or reversed interval" );
        if( i )
            require( records[ i - 1 ].validity.endTai && interval.startTai && *records[ i - 1 ].validity.endTai == *interval.startTai,
                     "history contains gaps or overlaps" );
    }
}
Eigen::Matrix3d localBasis( double latitude, double longitude )
{
    return reference_frames::getEnuLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion( longitude, latitude )
            .toRotationMatrix( );
}
void checkGeometry( const StationReferencePointGeometry& geometry )
{
    require( geometry.vertical == StationLocalVertical::explicitFixedBasis || geometry.vertical == StationLocalVertical::geocentricRadial ||
                     geometry.vertical == StationLocalVertical::geodeticEllipsoidNormal,
             "unknown local vertical convention" );
    const auto& matrix = geometry.enuToTerrestrial;
    require( !geometry.frame.empty( ) && matrix.allFinite( ) &&
                     ( matrix.transpose( ) * matrix - Eigen::Matrix3d::Identity( ) ).norm( ) < 1e-12 &&
                     std::abs( matrix.determinant( ) - 1 ) < 1e-12,
             "explicit geometry must be a proper fixed ENU rotation in a named frame" );
    if( geometry.vertical != StationLocalVertical::explicitFixedBasis )
    {
        const auto expected = geometry.vertical == StationLocalVertical::geocentricRadial
                ? createGeocentricStationReferencePointGeometry( geometry.anchorMetres, geometry.frame )
                : createGeodeticStationReferencePointGeometry(
                          geometry.anchorMetres,
                          geometry.frame,
                          geometry.equatorialRadiusMetres.value_or( std::numeric_limits< double >::quiet_NaN( ) ),
                          geometry.flattening.value_or( std::numeric_limits< double >::quiet_NaN( ) ) );
        require( ( matrix - expected.enuToTerrestrial ).norm( ) < 1e-12, "declared vertical/shape disagrees with supplied local basis" );
    }
}
std::optional< Time > boundary( const SinexEpoch& epoch,
                                bool start,
                                const SinexStationHistory& product,
                                const std::optional< SinexEquipmentHeaderTimeContext >& context )
{
    if( epoch.meaning == ( start ? SinexEpochMeaning::openStart : SinexEpochMeaning::openEnd ) ) return std::nullopt;
    const SinexEpoch* resolved = &epoch;
    if( epoch.meaning == ( start ? SinexEpochMeaning::headerStartDefault : SinexEpochMeaning::headerEndDefault ) )
    {
        resolved = start ? &product.headerStart : &product.headerEnd;
        if( resolved->meaning != SinexEpochMeaning::specified && context )
        {
            const auto& supplied = start ? context->startTai : context->endTai;
            require( !context->source.empty( ) && supplied.has_value( ), "missing explicit bounded header-time context" );
            return supplied;
        }
    }
    require( resolved->meaning == SinexEpochMeaning::specified && resolved->civilUtc.has_value( ),
             product.source + ": unresolved equipment header/default endpoint" );
    const Time tai = convertSinexCivilUtcToTai( *resolved->civilUtc );
    if( resolved != &epoch && context )
    {
        const auto& supplied = start ? context->startTai : context->endTai;
        require( !supplied || *supplied == tai, "header-time context contradicts a specified source header" );
    }
    return tai;
}
std::shared_ptr< StationReferencePointTrajectory > pointSource( const std::shared_ptr< StationHistoryTrajectory >& source )
{
    if( const auto point = std::dynamic_pointer_cast< StationReferencePointTrajectory >( source ) ) return point;
    if( const auto transformed = std::dynamic_pointer_cast< TransformedStationHistoryTrajectory >( source ) )
        return pointSource( transformed->getSourceTrajectory( ) );
    return nullptr;
}
}  // namespace

StationReferencePointGeometry createGeocentricStationReferencePointGeometry( const Eigen::Vector3d& anchor, const std::string& frame )
{
    require( anchor.allFinite( ) && anchor.norm( ) > 0 && anchor.head< 2 >( ).norm( ) > 0 && !frame.empty( ),
             "a finite nonpolar anchor and explicit frame are required" );
    StationReferencePointGeometry result;
    result.frame = frame;
    result.vertical = StationLocalVertical::geocentricRadial;
    result.anchorMetres = anchor;
    result.enuToTerrestrial = localBasis( std::atan2( anchor.z( ), anchor.head< 2 >( ).norm( ) ), std::atan2( anchor.y( ), anchor.x( ) ) );
    return result;
}
StationReferencePointGeometry createGeodeticStationReferencePointGeometry( const Eigen::Vector3d& anchor,
                                                                           const std::string& frame,
                                                                           double radius,
                                                                           double flattening )
{
    require( std::isfinite( radius ) && radius > 0 && std::isfinite( flattening ) && flattening >= 0 && flattening < 1,
             "positive equatorial radius and physical flattening are required" );
    auto result = createGeocentricStationReferencePointGeometry( anchor, frame );
    const Eigen::Vector3d coordinates = coordinate_conversions::convertCartesianToGeodeticCoordinates( anchor, radius, flattening, 1e-8 );
    result.vertical = StationLocalVertical::geodeticEllipsoidNormal;
    result.equatorialRadiusMetres = radius;
    result.flattening = flattening;
    result.enuToTerrestrial = localBasis( coordinates.y( ), coordinates.z( ) );
    return result;
}

StationReferencePointOffsetHistory::StationReferencePointOffsetHistory( std::vector< StationReferencePointOffset > offsets,
                                                                        std::optional< StationReferencePointGeometry > geometry,
                                                                        StationReferencePointRatePolicy policy ):
    offsets_( std::move( offsets ) ), geometry_( std::move( geometry ) ), ratePolicy_( policy )
{
    sortAndCheck( offsets_ );
    require( policy == StationReferencePointRatePolicy::rejectAtSwitch || policy == StationReferencePointRatePolicy::rightHandAtSwitch,
             "unknown switch-rate policy" );
    if( geometry_ ) checkGeometry( *geometry_ );
    for( const auto& offset : offsets_ )
    {
        for( const auto& id : { offset.from, offset.to } )
            require( !id.site.siteCode.empty( ) && !id.site.pointCode.empty( ) && !id.domesId.empty( ) && !id.name.empty( ),
                     "incomplete reference-point identity" );
        require( sameIdentity( offset.from, offsets_.front( ).from ) && sameIdentity( offset.to, offsets_.front( ).to ) &&
                         sameSite( offset.from.site, offset.to.site ) && offset.from.domesId == offset.to.domesId,
                 "offset history mixes physical point identities" );
        require( !offset.frame.empty( ) && offset.frame == offsets_.front( ).frame && offset.from.name != offset.to.name,
                 "explicit common frame and distinct point names required" );
        require( ( offset.from.kind == StationReferencePointKind::marker &&
                   offset.to.kind == StationReferencePointKind::antennaReferencePoint ) ||
                         ( offset.from.kind == StationReferencePointKind::antennaReferencePoint &&
                           offset.to.kind == StationReferencePointKind::calibratedInstrumentPoint ),
                 "only explicit marker-to-ARP or ARP-to-calibrated-point vectors are supported" );
        require( !offset.equipmentType.empty( ) && !offset.equipmentSerialNumber.empty( ) &&
                         ( !offset.provenance.source.empty( ) || !offset.source.empty( ) ) && offset.componentsMetres.allFinite( ) &&
                         !offset.directionDependent,
                 "incomplete equipment/source, nonfinite vector or direction-dependent correction" );
        require( offset.basis == StationReferencePointOffsetBasis::xyz ||
                         ( offset.basis == StationReferencePointOffsetBasis::une && geometry_ ),
                 "UNE offsets require explicit fixed local geometry" );
        require( !geometry_ || geometry_->frame == offset.frame, "offset component frame disagrees with its local geometry" );
        require( !offset.stateFunction || offset.positionFunction, "time-dependent state callback also needs a position callback" );
        if( offset.to.kind == StationReferencePointKind::calibratedInstrumentPoint )
            require( !offset.calibrationId.empty( ) && !offset.signal.empty( ),
                     "calibrated reference point requires explicit calibration and signal" );
    }
}
const StationReferencePointOffset& StationReferencePointOffsetHistory::getActiveOffset( const Time& tai ) const
{
    require( std::isfinite( static_cast< long double >( tai ) ), "nonfinite physical query epoch" );
    for( const auto& offset : offsets_ )
        if( contains( offset.validity, tai ) ) return offset;
    throw std::runtime_error( "Station reference point requested outside offset/equipment validity" );
}
Eigen::Vector3d StationReferencePointOffsetHistory::toTerrestrial( const Eigen::Vector3d& components,
                                                                   StationReferencePointOffsetBasis basis ) const
{
    require( components.allFinite( ), "nonfinite reference-point vector/rate" );
    if( basis == StationReferencePointOffsetBasis::xyz ) return components;
    return geometry_->enuToTerrestrial * Eigen::Vector3d( components.z( ), components.y( ), components.x( ) );
}
Eigen::Vector3d StationReferencePointOffsetHistory::getPositionCorrection( const StationHistoryTime& time ) const
{
    const auto& offset = getActiveOffset( time.tai );
    return toTerrestrial( offset.positionFunction ? offset.positionFunction( time ) : offset.componentsMetres, offset.basis );
}
Eigen::Vector6d StationReferencePointOffsetHistory::getStateCorrection( const StationHistoryTime& time ) const
{
    const auto& offset = getActiveOffset( time.tai );
    if( ratePolicy_ == StationReferencePointRatePolicy::rejectAtSwitch && &offset != &offsets_.front( ) && offset.validity.startTai &&
        time.tai == *offset.validity.startTai )
        throw std::runtime_error(
                "Reference-point state derivative is unavailable at an equipment/offset switch; choose explicit right-hand policy" );
    Eigen::Vector6d state = Eigen::Vector6d::Zero( );
    if( offset.positionFunction )
    {
        require( static_cast< bool >( offset.stateFunction ), "time-dependent reference-point rate is unavailable" );
        state = offset.stateFunction( time );
        const Eigen::Vector3d position = offset.positionFunction( time );
        require( state.allFinite( ) && position.allFinite( ) &&
                         ( state.head< 3 >( ) - position ).norm( ) <=
                                 32 * std::numeric_limits< double >::epsilon( ) * std::max( 1.0, position.norm( ) ),
                 "position/state calibration callbacks disagree" );
    }
    else
        state.head< 3 >( ) = offset.componentsMetres;
    state.head< 3 >( ) = toTerrestrial( state.head< 3 >( ), offset.basis );
    state.tail< 3 >( ) = toTerrestrial( state.tail< 3 >( ), offset.basis );
    return state;
}

//! Read the active reference-point identity from the trajectory's output metadata.
StationReferencePointIdentity getStationReferencePointIdentity( const std::shared_ptr< StationHistoryTrajectory >& source )
{
    require( source != nullptr, "a source trajectory is required" );

    // Output metadata preserves point identity through frame transformations.
    const auto& meta = source->getMetadata( );

    return { meta.site, meta.domesId, meta.referencePointKind, meta.referencePointName };
}
//! Retain the source trajectory and bind one validated reference-point offset history.
StationReferencePointTrajectory::StationReferencePointTrajectory( std::shared_ptr< StationHistoryTrajectory > source,
                                                                  std::shared_ptr< StationReferencePointOffsetHistory > offsets ):
    StationHistoryTrajectory( std::move( source ) ),
    offsets_( std::move( offsets ) ), outputMetadata_( source_->getMetadata( ) )
{
    require( offsets_ != nullptr, "an explicit reference-point offset history is required" );
    require( sameIdentity( getStationReferencePointIdentity( source_ ), offsets_->getOffsets( ).front( ).from ),
             "source point identity disagrees with offset direction; duplicate offsets are not permitted" );
    require( !offsets_->getGeometry( ) || offsets_->getGeometry( )->frame == outputMetadata_.frame,
             "offset local basis and trajectory frames differ" );
    require( offsets_->getOffsets( ).front( ).frame == outputMetadata_.frame, "offset and source trajectory frames differ" );
    outputMetadata_.referencePointKind = getReferencePointIdentity( ).kind;
    outputMetadata_.referencePointName = getReferencePointIdentity( ).name;
    for( const auto& record : offsets_->getOffsets( ) )
    {
        if( record.to.kind == StationReferencePointKind::calibratedInstrumentPoint )
        {
            // A calibration must remain inside the parent ARP equipment coverage and use its equipment identity.
            const auto parent = pointSource( source_ );

            require( parent != nullptr, "calibration requires an explicitly bound ARP equipment history" );
            const auto& parentRecords = parent->getOffsetHistory( )->getOffsets( );

            require( ( !parentRecords.front( ).validity.startTai ||
                       ( record.validity.startTai && *record.validity.startTai >= *parentRecords.front( ).validity.startTai ) ) &&
                             ( !parentRecords.back( ).validity.endTai ||
                               ( record.validity.endTai && *record.validity.endTai <= *parentRecords.back( ).validity.endTai ) ),
                     "calibration extends beyond equipment validity" );
            for( const auto& equipment : parentRecords )
            {
                if( overlap( record.validity, equipment.validity ) )
                {
                    require( record.equipmentType == equipment.equipmentType &&
                                     record.equipmentSerialNumber == equipment.equipmentSerialNumber,
                             "calibration contradicts the active equipment identity" );
                }
            }
        }
        for( const auto& source : { record.source,
                                    record.provenance.source,
                                    record.equipmentProvenance.source,
                                    record.headerTimeContext ? record.headerTimeContext->source : std::string( ) } )
        {
            if( !source.empty( ) &&
                std::find( outputMetadata_.sources.begin( ), outputMetadata_.sources.end( ), source ) == outputMetadata_.sources.end( ) )
            {
                outputMetadata_.sources.push_back( source );
            }
        }
    }
    outputMetadata_.signals.push_back( getReferencePointIdentity( ).kind == StationReferencePointKind::antennaReferencePoint
                                               ? "marker-to-antenna reference vector"
                                               : "calibrated instrument reference vector" );
}
//! Add the active offset in the declared terrestrial frame to the source position.
Eigen::Vector3d StationReferencePointTrajectory::getPosition( const StationHistoryTime& time ) const
{
    // Evaluate source motion first, preserving the explicitly selected composition order.
    Eigen::Vector3d position = source_->getPosition( time );

    position += offsets_->getPositionCorrection( time );
    return position;
}
//! Add the active offset and its TDB velocity to the complete source state.
Eigen::Vector6d StationReferencePointTrajectory::getState( const StationHistoryTime& time ) const
{
    // Source and offset velocities are derivatives with respect to the same TDB seconds.
    Eigen::Vector6d state = source_->getState( time );

    state += offsets_->getStateCorrection( time );
    return state;
}
std::shared_ptr< StationReferencePointTrajectory > createStationReferencePointTrajectory(
        std::shared_ptr< StationHistoryTrajectory > source,
        std::shared_ptr< StationReferencePointOffsetHistory > offsets )
{
    return std::make_shared< StationReferencePointTrajectory >( std::move( source ), std::move( offsets ) );
}

std::shared_ptr< StationReferencePointTrajectory > createSinexAntennaReferencePointTrajectory(
        std::shared_ptr< StationHistoryTrajectory > marker,
        const SinexStationHistory& equipment,
        std::string name,
        std::optional< StationReferencePointGeometry > geometry,
        std::optional< SinexEquipmentHeaderTimeContext > context,
        StationReferencePointRatePolicy policy,
        const std::optional< std::string >& domesSelector,
        const std::optional< SinexSiteKey >& equipmentSite )
{
    validateSinexStationHistory( equipment );
    const auto from = getStationReferencePointIdentity( marker );
    require( from.kind == StationReferencePointKind::marker && !name.empty( ), "an explicit marker and ARP name are required" );
    const auto& meta = marker->getMetadata( );
    const auto key = equipmentSite.value_or( from.site );
    require( equipment.frame.empty( ) || equipment.frame == meta.frame, "equipment/marker coordinate frames differ" );
    if( context )
    {
        require( !context->source.empty( ), "explicit header-time context requires provenance" );
        for( const auto& value : { context->startTai, context->endTai } )
            require( !value || std::isfinite( static_cast< long double >( *value ) ), "nonfinite equipment header context" );
        require( !context->startTai || !context->endTai || *context->startTai < *context->endTai, "reversed equipment header context" );
    }
    std::set< std::string > domes;
    const auto sites = equipment.sites.equal_range( key );
    for( auto it = sites.first; it != sites.second; ++it ) domes.insert( it->second.domesId );
    require( domes.count( from.domesId ) && ( domesSelector ? *domesSelector == from.domesId : domes.size( ) == 1 ),
             "missing or ambiguous equipment DOMES association; explicit matching selector required" );
    const StationReferencePointIdentity to{ from.site, from.domesId, StationReferencePointKind::antennaReferencePoint, std::move( name ) };
    struct Eccentricity {
        SinexEccentricityRecord source;
        StationHistoryInterval validity;
    };
    struct Antenna {
        SinexEquipmentRecord source;
        StationHistoryInterval validity;
    };
    std::vector< Eccentricity > eccentricities;
    std::vector< Antenna > antennas;
    for( const auto& record : equipment.eccentricities )
        if( sameSite( record.key.site, key ) )
        {
            require( record.observationTechnique == "P", "this factory interprets GNSS marker-to-ARP eccentricity only" );
            if( record.basis == SinexEccentricityBasis::xyz )
                require( equipment.frame == meta.frame, "XYZ eccentricity requires an explicitly matching coordinate frame" );
            eccentricities.push_back(
                    { record, { boundary( record.start, true, equipment, context ), boundary( record.end, false, equipment, context ) } } );
        }
    for( const auto& record : equipment.equipmentHistories )
        if( sameSite( record.key.site, key ) && record.kind == SinexEquipmentKind::antenna )
        {
            require( record.observationTechnique == "P", "GNSS antenna equipment association required" );
            antennas.push_back(
                    { record, { boundary( record.start, true, equipment, context ), boundary( record.end, false, equipment, context ) } } );
        }
    sortAndCheck( eccentricities );
    sortAndCheck( antennas );
    require( eccentricities.front( ).validity.startTai == antennas.front( ).validity.startTai &&
                     eccentricities.back( ).validity.endTai == antennas.back( ).validity.endTai,
             "antenna/eccentricity coverage differs" );
    std::vector< StationReferencePointOffset > offsets;
    for( const auto& eccentricity : eccentricities )
        for( const auto& antenna : antennas )
            if( overlap( eccentricity.validity, antenna.validity ) )
            {
                require( eccentricity.source.key.solutionId == antenna.source.key.solutionId,
                         "equipment/eccentricity source solution identifiers differ" );
                StationReferencePointOffset record;
                record.from = from;
                record.to = to;
                record.frame = meta.frame;
                record.validity = intersection( eccentricity.validity, antenna.validity );
                record.basis = eccentricity.source.basis == SinexEccentricityBasis::une ? StationReferencePointOffsetBasis::une
                                                                                        : StationReferencePointOffsetBasis::xyz;
                for( int i = 0; i < 3; ++i ) record.componentsMetres[ i ] = eccentricity.source.componentsMetres[ i ];
                record.equipmentType = antenna.source.equipmentType;
                record.equipmentSerialNumber = antenna.source.serialNumber;
                record.provenance = eccentricity.source.provenance;
                record.equipmentProvenance = antenna.source.provenance;
                record.headerTimeContext = context;
                offsets.push_back( record );
            }
    auto result = createStationReferencePointTrajectory(
            marker, std::make_shared< StationReferencePointOffsetHistory >( offsets, geometry, policy ) );
    return result;
}
}  // namespace ground_stations
}  // namespace tudat

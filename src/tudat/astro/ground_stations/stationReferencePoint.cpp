/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/astro/ground_stations/stationReferencePoint.h"
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
void requireValidReferencePoint( bool condition, const std::string& message )
{
    if( !condition )
    {
        throw std::runtime_error( "Station reference point: " + message );
    }
}
bool sameIdentity( const StationReferencePointIdentity& a, const StationReferencePointIdentity& b )
{
    return a.site == b.site && a.domesId == b.domesId && a.kind == b.kind && a.name == b.name;
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
    requireValidReferencePoint( !records.empty( ), "empty offset/equipment history" );
    std::sort( records.begin( ), records.end( ), []( const auto& a, const auto& b ) {
        return !a.validity.startTai ? b.validity.startTai.has_value( )
                                    : ( b.validity.startTai && *a.validity.startTai < *b.validity.startTai );
    } );
    for( std::size_t i = 0; i < records.size( ); ++i )
    {
        const auto& interval = records[ i ].validity;
        for( const auto& boundary : { interval.startTai, interval.endTai } )
        {
            requireValidReferencePoint( !boundary || std::isfinite( static_cast< long double >( *boundary ) ),
                                        "nonfinite validity endpoint" );
        }
        requireValidReferencePoint( !interval.startTai || !interval.endTai || *interval.startTai < *interval.endTai,
                                    "empty or reversed interval" );
        if( i )
        {
            requireValidReferencePoint(
                    records[ i - 1 ].validity.endTai && interval.startTai && *records[ i - 1 ].validity.endTai == *interval.startTai,
                    "history contains gaps or overlaps" );
        }
    }
}
Eigen::Matrix3d localBasis( double latitude, double longitude )
{
    return reference_frames::getEnuLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion( longitude, latitude )
            .toRotationMatrix( );
}
void checkGeometry( const StationReferencePointGeometry& geometry )
{
    requireValidReferencePoint( geometry.vertical == StationLocalVertical::explicitFixedBasis ||
                                        geometry.vertical == StationLocalVertical::geocentricRadial ||
                                        geometry.vertical == StationLocalVertical::geodeticEllipsoidNormal,
                                "unknown local vertical convention" );
    const auto& matrix = geometry.enuToTerrestrial;
    requireValidReferencePoint( !geometry.frame.empty( ) && matrix.allFinite( ) &&
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
        requireValidReferencePoint( ( matrix - expected.enuToTerrestrial ).norm( ) < 1e-12,
                                    "declared vertical/shape disagrees with supplied local basis" );
    }
}
std::optional< Time > boundary( const SinexEpoch& epoch,
                                bool start,
                                const SinexStationHistory& product,
                                const std::optional< SinexEquipmentHeaderTimeContext >& context )
{
    if( epoch.meaning == ( start ? SinexEpochMeaning::openStart : SinexEpochMeaning::openEnd ) )
    {
        return std::nullopt;
    }
    const SinexEpoch* resolved = &epoch;
    if( epoch.meaning == ( start ? SinexEpochMeaning::headerStartDefault : SinexEpochMeaning::headerEndDefault ) )
    {
        resolved = start ? &product.headerStart : &product.headerEnd;
        if( resolved->meaning != SinexEpochMeaning::specified && context )
        {
            const auto& supplied = start ? context->startTai : context->endTai;
            requireValidReferencePoint( !context->source.empty( ) && supplied.has_value( ),
                                        "missing explicit bounded header-time context" );
            return supplied;
        }
    }
    requireValidReferencePoint( resolved->meaning == SinexEpochMeaning::specified && resolved->civilUtc.has_value( ),
                                product.source + ": unresolved equipment header/default endpoint" );
    const Time tai = convertSinexCivilUtcToTai( *resolved->civilUtc );
    if( resolved != &epoch && context )
    {
        const auto& supplied = start ? context->startTai : context->endTai;
        requireValidReferencePoint( !supplied || *supplied == tai, "header-time context contradicts a specified source header" );
    }
    return tai;
}
}  // namespace

StationReferencePointGeometry createGeocentricStationReferencePointGeometry( const Eigen::Vector3d& anchor, const std::string& frame )
{
    requireValidReferencePoint( anchor.allFinite( ) && anchor.norm( ) > 0 && anchor.head< 2 >( ).norm( ) > 0 && !frame.empty( ),
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
    requireValidReferencePoint( std::isfinite( radius ) && radius > 0 && std::isfinite( flattening ) && flattening >= 0 && flattening < 1,
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
                                                                        std::optional< StationReferencePointGeometry > geometry ):
    offsets_( std::move( offsets ) ), geometry_( std::move( geometry ) )
{
    sortAndCheck( offsets_ );
    if( geometry_ )
    {
        checkGeometry( *geometry_ );
    }
    for( const auto& offset : offsets_ )
    {
        for( const auto& id : { offset.from, offset.to } )
        {
            requireValidReferencePoint(
                    !id.site.siteCode.empty( ) && !id.site.pointCode.empty( ) && !id.domesId.empty( ) && !id.name.empty( ),
                    "incomplete reference-point identity" );
        }
        requireValidReferencePoint( sameIdentity( offset.from, offsets_.front( ).from ) &&
                                            sameIdentity( offset.to, offsets_.front( ).to ) && offset.from.site == offset.to.site &&
                                            offset.from.domesId == offset.to.domesId,
                                    "offset history mixes physical point identities" );
        requireValidReferencePoint( !offset.frame.empty( ) && offset.frame == offsets_.front( ).frame && offset.from.name != offset.to.name,
                                    "explicit common frame and distinct point names required" );
        requireValidReferencePoint( ( offset.from.kind == StationReferencePointKind::marker &&
                                      offset.to.kind == StationReferencePointKind::antennaReferencePoint ) ||
                                            ( offset.from.kind == StationReferencePointKind::antennaReferencePoint &&
                                              offset.to.kind == StationReferencePointKind::calibratedInstrumentPoint ),
                                    "only explicit marker-to-ARP or ARP-to-calibrated-point vectors are supported" );
        requireValidReferencePoint( !offset.equipmentType.empty( ) && !offset.equipmentSerialNumber.empty( ) &&
                                            ( !offset.provenance.source.empty( ) || !offset.source.empty( ) ) &&
                                            offset.componentsMetres.allFinite( ) && !offset.directionDependent,
                                    "incomplete equipment/source, nonfinite vector or direction-dependent correction" );
        requireValidReferencePoint( offset.basis == StationReferencePointOffsetBasis::xyz ||
                                            ( offset.basis == StationReferencePointOffsetBasis::une && geometry_ ),
                                    "UNE offsets require explicit fixed local geometry" );
        requireValidReferencePoint( !geometry_ || geometry_->frame == offset.frame,
                                    "offset component frame disagrees with its local geometry" );
        requireValidReferencePoint( !offset.stateFunction || offset.positionFunction,
                                    "time-dependent state callback also needs a position callback" );
        if( offset.to.kind == StationReferencePointKind::calibratedInstrumentPoint )
        {
            requireValidReferencePoint( !offset.calibrationId.empty( ) && !offset.signal.empty( ),
                                        "calibrated reference point requires explicit calibration and signal" );
        }
    }
}
const StationReferencePointOffset& StationReferencePointOffsetHistory::getActiveOffset( const Time& tai ) const
{
    requireValidReferencePoint( std::isfinite( static_cast< long double >( tai ) ), "nonfinite physical query epoch" );
    for( const auto& offset : offsets_ )
    {
        if( contains( offset.validity, tai ) )
        {
            return offset;
        }
    }
    throw std::runtime_error( "Station reference point requested outside offset/equipment validity" );
}
Eigen::Vector3d StationReferencePointOffsetHistory::toTerrestrial( const Eigen::Vector3d& components,
                                                                   StationReferencePointOffsetBasis basis ) const
{
    requireValidReferencePoint( components.allFinite( ), "nonfinite reference-point vector/rate" );
    if( basis == StationReferencePointOffsetBasis::xyz )
    {
        return components;
    }
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
    Eigen::Vector6d state = Eigen::Vector6d::Zero( );
    if( offset.positionFunction )
    {
        requireValidReferencePoint( static_cast< bool >( offset.stateFunction ), "time-dependent reference-point rate is unavailable" );
        state = offset.stateFunction( time );
        const Eigen::Vector3d position = offset.positionFunction( time );
        requireValidReferencePoint( state.allFinite( ) && position.allFinite( ) &&
                                            ( state.head< 3 >( ) - position ).norm( ) <=
                                                    32 * std::numeric_limits< double >::epsilon( ) * std::max( 1.0, position.norm( ) ),
                                    "position/state calibration callbacks disagree" );
    }
    else
    {
        state.head< 3 >( ) = offset.componentsMetres;
    }
    state.head< 3 >( ) = toTerrestrial( state.head< 3 >( ), offset.basis );
    state.tail< 3 >( ) = toTerrestrial( state.tail< 3 >( ), offset.basis );
    return state;
}

//! Identity of the point described by a station coordinate history.
StationReferencePointIdentity getStationReferencePointIdentity( const StationHistoryMetadata& metadata )
{
    return { metadata.site, metadata.domesId, metadata.referencePointKind, metadata.referencePointName };
}

//! Check that the offset histories form a chain starting at the point of the history, in the frame of the history.
void validateStationReferencePointOffsets( const StationHistoryMetadata& history,
                                           const std::vector< std::shared_ptr< StationReferencePointOffsetHistory > >& offsets )
{
    StationReferencePointIdentity previousPoint = getStationReferencePointIdentity( history );
    std::shared_ptr< StationReferencePointOffsetHistory > previous;
    for( const auto& current : offsets )
    {
        requireValidReferencePoint( current != nullptr, "an explicit reference-point offset history is required" );
        const auto& records = current->getOffsets( );
        requireValidReferencePoint( sameIdentity( previousPoint, records.front( ).from ),
                                    "offsets do not start at the point described by the station history (or the previous offsets)" );
        requireValidReferencePoint(
                records.front( ).frame == history.frame && ( !current->getGeometry( ) || current->getGeometry( )->frame == history.frame ),
                "offset and station history frames differ" );
        for( const auto& record : records )
        {
            if( record.to.kind != StationReferencePointKind::calibratedInstrumentPoint )
            {
                continue;
            }

            // A calibration must remain inside the coverage of the antenna offsets before it, and use the same antenna.
            requireValidReferencePoint( previous != nullptr, "calibration requires an explicitly bound ARP equipment history" );
            const auto& parentRecords = previous->getOffsets( );
            requireValidReferencePoint(
                    ( !parentRecords.front( ).validity.startTai ||
                      ( record.validity.startTai && *record.validity.startTai >= *parentRecords.front( ).validity.startTai ) ) &&
                            ( !parentRecords.back( ).validity.endTai ||
                              ( record.validity.endTai && *record.validity.endTai <= *parentRecords.back( ).validity.endTai ) ),
                    "calibration extends beyond equipment validity" );
            for( const auto& equipment : parentRecords )
            {
                if( overlap( record.validity, equipment.validity ) )
                {
                    requireValidReferencePoint( record.equipmentType == equipment.equipmentType &&
                                                        record.equipmentSerialNumber == equipment.equipmentSerialNumber,
                                                "calibration contradicts the active equipment identity" );
                }
            }
        }
        previousPoint = records.front( ).to;
        previous = current;
    }
}

std::shared_ptr< StationReferencePointOffsetHistory > createSinexAntennaReferencePointOffsets(
        const StationHistoryMetadata& meta,
        const SinexStationHistory& equipment,
        std::string name,
        std::optional< StationReferencePointGeometry > geometry,
        std::optional< SinexEquipmentHeaderTimeContext > context,
        const std::optional< std::string >& domesSelector,
        const std::optional< SinexSiteKey >& equipmentSite )
{
    validateSinexStationHistory( equipment );
    const auto from = getStationReferencePointIdentity( meta );
    requireValidReferencePoint( from.kind == StationReferencePointKind::marker && !name.empty( ),
                                "an explicit marker and ARP name are required" );
    const auto key = equipmentSite.value_or( from.site );
    requireValidReferencePoint( equipment.frame.empty( ) || equipment.frame == meta.frame, "equipment/marker coordinate frames differ" );
    if( context )
    {
        requireValidReferencePoint( !context->source.empty( ), "explicit header-time context requires provenance" );
        for( const auto& value : { context->startTai, context->endTai } )
        {
            requireValidReferencePoint( !value || std::isfinite( static_cast< long double >( *value ) ),
                                        "nonfinite equipment header context" );
        }
        requireValidReferencePoint( !context->startTai || !context->endTai || *context->startTai < *context->endTai,
                                    "reversed equipment header context" );
    }
    std::set< std::string > domes;
    const auto sites = equipment.sites.equal_range( key );
    for( auto it = sites.first; it != sites.second; ++it )
    {
        domes.insert( it->second.domesId );
    }
    requireValidReferencePoint( domes.count( from.domesId ) && ( domesSelector ? *domesSelector == from.domesId : domes.size( ) == 1 ),
                                "missing or ambiguous equipment DOMES association; explicit matching selector required" );
    const StationReferencePointIdentity to{ from.site, from.domesId, StationReferencePointKind::antennaReferencePoint, std::move( name ) };

    //! Eccentricity record of the station, with its validity interval in TAI.
    struct Eccentricity {
        //! Record as read, and its validity interval.
        SinexEccentricityRecord source;
        StationHistoryInterval validity;
    };

    //! Antenna record of the station, with its installation interval in TAI.
    struct Antenna {
        //! Record as read, and its installation interval.
        SinexEquipmentRecord source;
        StationHistoryInterval validity;
    };
    std::vector< Eccentricity > eccentricities;
    std::vector< Antenna > antennas;
    for( const auto& record : equipment.eccentricities )
    {
        if( record.key.site == key )
        {
            requireValidReferencePoint( record.observationTechnique == "P",
                                        "this factory interprets GNSS marker-to-ARP eccentricity only" );
            if( record.basis == SinexEccentricityBasis::xyz )
            {
                requireValidReferencePoint( equipment.frame == meta.frame,
                                            "XYZ eccentricity requires an explicitly matching coordinate frame" );
            }
            eccentricities.push_back(
                    { record, { boundary( record.start, true, equipment, context ), boundary( record.end, false, equipment, context ) } } );
        }
    }
    for( const auto& record : equipment.equipmentHistories )
    {
        if( record.key.site == key && record.kind == SinexEquipmentKind::antenna )
        {
            requireValidReferencePoint( record.observationTechnique == "P", "GNSS antenna equipment association required" );
            antennas.push_back(
                    { record, { boundary( record.start, true, equipment, context ), boundary( record.end, false, equipment, context ) } } );
        }
    }
    sortAndCheck( eccentricities );
    sortAndCheck( antennas );
    requireValidReferencePoint( eccentricities.front( ).validity.startTai == antennas.front( ).validity.startTai &&
                                        eccentricities.back( ).validity.endTai == antennas.back( ).validity.endTai,
                                "antenna/eccentricity coverage differs" );
    std::vector< StationReferencePointOffset > offsets;
    for( const auto& eccentricity : eccentricities )
    {
        for( const auto& antenna : antennas )
        {
            if( overlap( eccentricity.validity, antenna.validity ) )
            {
                requireValidReferencePoint( eccentricity.source.key.solutionId == antenna.source.key.solutionId,
                                            "equipment/eccentricity source solution identifiers differ" );
                StationReferencePointOffset record;
                record.from = from;
                record.to = to;
                record.frame = meta.frame;
                record.validity = intersection( eccentricity.validity, antenna.validity );
                record.basis = eccentricity.source.basis == SinexEccentricityBasis::une ? StationReferencePointOffsetBasis::une
                                                                                        : StationReferencePointOffsetBasis::xyz;
                for( int i = 0; i < 3; ++i )
                {
                    record.componentsMetres[ i ] = eccentricity.source.componentsMetres[ i ];
                }
                record.equipmentType = antenna.source.equipmentType;
                record.equipmentSerialNumber = antenna.source.serialNumber;
                record.provenance = eccentricity.source.provenance;
                record.equipmentProvenance = antenna.source.provenance;
                record.headerTimeContext = context;
                offsets.push_back( record );
            }
        }
    }
    return std::make_shared< StationReferencePointOffsetHistory >( offsets, geometry );
}
}  // namespace ground_stations
}  // namespace tudat

/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/astro/ground_stations/stationHistoryTrajectory.h"
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

namespace tudat
{
namespace ground_stations
{
namespace
{
using namespace input_output;
constexpr double daysPerJulianYear = 365.25;
constexpr double secondsPerJulianYear = 31557600.0;
bool sameSite( const SinexSiteKey& a, const SinexSiteKey& b )
{
    return !( a < b ) && !( b < a );
}
bool sameKey( const SinexSolutionKey& a, const SinexSolutionKey& b )
{
    return !( a < b ) && !( b < a );
}
void require( bool condition, const std::string& message )
{
    if( !condition ) throw std::runtime_error( "Station history: " + message );
}
void validateTime( const StationHistoryTime& time, const std::string& convention )
{
    require( std::isfinite( static_cast< long double >( time.tai ) ) && std::isfinite( time.mjdFraction ) && time.mjdFraction >= 0.0L &&
                     time.mjdFraction < 1.0L && std::isfinite( time.yearFraction ) && time.yearFraction >= 0.0L && time.yearFraction < 1.0L,
             "invalid split source coordinates" );
    require( !convention.empty( ) && time.convention == convention, "source time conventions differ" );
    for( const auto& rate : { time.mjdRate, time.decimalYearRate } )
        require( !rate || ( std::isfinite( *rate ) && *rate > 0.0 ), "source coordinate rate must be finite and positive" );
}
double elapsedYears( const StationHistoryTime& time, const StationHistoryTime& reference )
{
    return static_cast< double >(
            ( static_cast< long double >( time.mjdDay ) - reference.mjdDay + time.mjdFraction - reference.mjdFraction ) /
            daysPerJulianYear );
}
std::string uniqueDomes( const SinexStationHistory& history,
                         const SinexSiteKey& site,
                         const std::optional< std::string >& selector = std::nullopt )
{
    const auto range = history.sites.equal_range( site );
    std::set< std::string > identifiers;
    for( auto it = range.first; it != range.second; ++it ) identifiers.insert( it->second.domesId );
    if( selector )
    {
        require( !selector->empty( ) && identifiers.count( *selector ), history.source + ": selected DOMES is absent" );
        return *selector;
    }
    require( identifiers.size( ) == 1 && !identifiers.begin( )->empty( ),
             history.source + ": station/point has missing or ambiguous DOMES association" );
    return *identifiers.begin( );
}
void checkMetadata( const SinexStationHistory& product, const StationHistoryMetadata& metadata )
{
    validateSinexStationHistory( product );
    require( !product.frame.empty( ) && !product.release.empty( ) && !product.origin.empty( ) && !product.tideSystem.empty( ),
             product.source + ": explicit frame, release, origin and tide-system metadata required" );
    require( product.frame == metadata.frame && product.release == metadata.release && product.origin == metadata.origin &&
                     product.tideSystem == metadata.tideSystem,
             product.source + ": mixed frame, release, origin or tide-system metadata" );
}
std::optional< Time > intervalEnd( const SinexEpoch& epoch, bool start, const SinexStationHistory& history )
{
    if( epoch.meaning == ( start ? SinexEpochMeaning::openStart : SinexEpochMeaning::openEnd ) ) return std::nullopt;
    const SinexEpoch* resolved = &epoch;
    if( epoch.meaning == ( start ? SinexEpochMeaning::headerStartDefault : SinexEpochMeaning::headerEndDefault ) )
        resolved = start ? &history.headerStart : &history.headerEnd;
    require( resolved->meaning == SinexEpochMeaning::specified && resolved->civilUtc.has_value( ),
             history.source + ": unresolved station validity endpoint" );
    return convertSinexCivilUtcToTai( *resolved->civilUtc );
}
StationHistoryTime sourceEpoch( const SinexEpoch& epoch )
{
    require( epoch.meaning == SinexEpochMeaning::specified && epoch.civilUtc.has_value( ), "unspecified model reference/event epoch" );
    return getIgnNominalUtcStationTime( *epoch.civilUtc );
}
Eigen::Vector3d componentVector( char component, const std::optional< Eigen::Matrix3d >& enuToTerrestrial )
{
    const auto xyz = std::string( "XYZ" ).find( component );
    if( xyz != std::string::npos ) return Eigen::Vector3d::Unit( xyz );
    const auto enu = std::string( "ENU" ).find( component );
    require( enu != std::string::npos && enuToTerrestrial.has_value( ), "ENU coefficients require an explicit constant local basis" );
    return enuToTerrestrial->col( enu );
}
double metres( const SinexParameterRecord& parameter )
{
    require( parameter.normalization == SinexNormalization::si && parameter.siUnit == "m" && parameter.valueSi.has_value( ),
             parameter.provenance.source + ": unsupported amplitude unit at parameter " + std::to_string( parameter.index ) );
    return *parameter.valueSi;
}
}  // namespace

StationHistoryTrajectory::StationHistoryTrajectory( std::vector< StationHistoryLinearSolution > solutions,
                                                    std::vector< StationPostSeismicTerm > postSeismicTerms,
                                                    std::vector< StationSeasonalTerm > seasonalTerms,
                                                    StationHistoryMetadata metadata ):
    solutions_( std::move( solutions ) ), postSeismicTerms_( std::move( postSeismicTerms ) ), seasonalTerms_( std::move( seasonalTerms ) ),
    metadata_( std::move( metadata ) )
{
    require( !solutions_.empty( ), "no linear position solutions supplied" );
    require( !metadata_.frame.empty( ) && !metadata_.release.empty( ), "explicit frame and release metadata required" );
    std::sort( solutions_.begin( ), solutions_.end( ), []( const auto& a, const auto& b ) {
        return !a.validity.startTai ? b.validity.startTai.has_value( )
                                    : ( b.validity.startTai && *a.validity.startTai < *b.validity.startTai );
    } );
    std::set< std::string > identifiers;
    for( std::size_t i = 0; i < solutions_.size( ); ++i )
    {
        const auto& solution = solutions_[ i ];
        const auto& interval = solution.validity;
        validateTime( solution.referenceTime, metadata_.timeConvention );
        require( solution.positionMetres.allFinite( ) && solution.velocityMetresPerJulianYear.allFinite( ), "nonfinite linear solution" );
        require( sameSite( solution.source.key.site, metadata_.site ) && !solution.source.key.solutionId.empty( ) &&
                         identifiers.insert( solution.source.key.solutionId ).second,
                 "inconsistent or repeated position solution identity" );
        for( const auto& boundary : { interval.startTai, interval.endTai } )
            require( !boundary || std::isfinite( static_cast< long double >( *boundary ) ), "nonfinite validity endpoint" );
        require( !interval.startTai || !interval.endTai || *interval.startTai < *interval.endTai, "empty or reversed validity interval" );
        if( i > 0 )
        {
            const auto& previous = solutions_[ i - 1 ].validity;
            require( previous.endTai && interval.startTai && *previous.endTai == *interval.startTai,
                     "position validity intervals overlap or contain a gap" );
        }
    }
    for( const auto& term : postSeismicTerms_ )
    {
        validateTime( term.eventTime, metadata_.timeConvention );
        require( term.amplitudeMetres.allFinite( ) && std::isfinite( term.relaxationJulianYears ) && term.relaxationJulianYears > 0.0,
                 "post-seismic relaxation must be positive and finite" );
        require( term.function == StationPostSeismicFunction::logarithmic || term.function == StationPostSeismicFunction::exponential,
                 "unsupported post-seismic function" );
    }
    std::set< std::pair< std::string, int > > seasonalKeys;
    for( const auto& term : seasonalTerms_ )
    {
        require( identifiers.count( term.solutionId ) && ( term.cyclesPerCalendarYear == 1 || term.cyclesPerCalendarYear == 2 ) &&
                         seasonalKeys.emplace( term.solutionId, term.cyclesPerCalendarYear ).second,
                 "unknown or repeated seasonal position-solution/frequency identity" );
        require( term.cosineMetres.allFinite( ) && term.sineMetres.allFinite( ), "nonfinite seasonal coefficient" );
    }
}

Eigen::Matrix< double, 6, 1 > StationHistoryTrajectory::evaluate( const StationHistoryTime& time, bool rates ) const
{
    validateTime( time, metadata_.timeConvention );
    const StationHistoryLinearSolution* selected = nullptr;
    for( const auto& solution : solutions_ )
        if( ( !solution.validity.startTai || time.tai >= *solution.validity.startTai ) &&
            ( !solution.validity.endTai || time.tai < *solution.validity.endTai ) )
            selected = &solution;
    require( selected != nullptr, "requested epoch outside physical station validity" );
    require( !rates || time.mjdRate.has_value( ), "MJD derivative with respect to TDB is unavailable" );
    const double yearsRate = rates ? *time.mjdRate / daysPerJulianYear : 0.0;
    Eigen::Matrix< double, 6, 1 > result = Eigen::Matrix< double, 6, 1 >::Zero( );
    result.head< 3 >( ) = selected->positionMetres + selected->velocityMetresPerJulianYear * elapsedYears( time, selected->referenceTime );
    if( rates ) result.tail< 3 >( ) = selected->velocityMetresPerJulianYear * yearsRate;
    for( const auto& term : postSeismicTerms_ )
    {
        if( time.tai < term.eventTime.tai ) continue;
        const double dt = elapsedYears( time, term.eventTime );
        require( dt >= 0.0, "physical and source event ordering disagree" );
        const double x = dt / term.relaxationJulianYears;
        if( term.function == StationPostSeismicFunction::logarithmic )
        {
            result.head< 3 >( ) += term.amplitudeMetres * std::log1p( x );
            if( rates ) result.tail< 3 >( ) += term.amplitudeMetres * ( yearsRate / ( term.relaxationJulianYears + dt ) );
        }
        else
        {
            result.head< 3 >( ) -= term.amplitudeMetres * std::expm1( -x );
            if( rates ) result.tail< 3 >( ) += term.amplitudeMetres * ( yearsRate / term.relaxationJulianYears * std::exp( -x ) );
        }
    }
    for( const auto& term : seasonalTerms_ )
    {
        if( term.solutionId != selected->source.key.solutionId ) continue;
        require( !rates || time.decimalYearRate.has_value( ), "calendar decimal-year derivative with respect to TDB is unavailable" );
        const double omega = 2.0 * std::acos( -1.0 ) * term.cyclesPerCalendarYear;
        const double angle = static_cast< double >( omega * time.yearFraction );
        result.head< 3 >( ) += term.cosineMetres * std::cos( angle ) + term.sineMetres * std::sin( angle );
        if( rates )
            result.tail< 3 >( ) +=
                    ( -term.cosineMetres * std::sin( angle ) + term.sineMetres * std::cos( angle ) ) * ( omega * *time.decimalYearRate );
    }
    require( result.allFinite( ), "nonfinite trajectory result" );
    return result;
}
Eigen::Vector3d StationHistoryTrajectory::getPosition( const StationHistoryTime& time ) const
{
    return evaluate( time, false ).head< 3 >( );
}
Eigen::Matrix< double, 6, 1 > StationHistoryTrajectory::getState( const StationHistoryTime& time ) const
{
    return evaluate( time, true );
}

std::shared_ptr< StationHistoryTrajectory > createIgnStationHistoryTrajectory( const SinexStationHistory& solutions,
                                                                               const SinexStationHistory& discontinuities,
                                                                               const SinexSiteKey& site,
                                                                               const SinexStationHistory* postSeismic,
                                                                               const SinexStationHistory* seasonal,
                                                                               const std::optional< Eigen::Matrix3d >& enuToTerrestrial,
                                                                               const std::optional< std::string >& domesSelector,
                                                                               double velocityConsistencyToleranceMetresPerJulianYear )
{
    require( std::isfinite( velocityConsistencyToleranceMetresPerJulianYear ) && velocityConsistencyToleranceMetresPerJulianYear >= 0.0,
             "invalid declared velocity consistency tolerance" );
    StationHistoryMetadata metadata;
    metadata.site = site;
    metadata.domesId = uniqueDomes( solutions, site, domesSelector );
    metadata.frame = solutions.frame;
    metadata.release = solutions.release;
    metadata.origin = solutions.origin;
    metadata.tideSystem = solutions.tideSystem;
    metadata.timeConvention = "IGN nominal UTC";
    metadata.signals = { "piecewise linear station coordinates" };
    for( const auto* product : { &solutions, &discontinuities, postSeismic, seasonal } )
    {
        if( !product ) continue;
        checkMetadata( *product, metadata );
        metadata.sources.push_back( product->source );
        if( product == postSeismic || product == seasonal )
            require( uniqueDomes( *product, site, domesSelector ) == metadata.domesId,
                     product->source + ": station monument association differs" );
    }
    if( enuToTerrestrial )
        require( enuToTerrestrial->allFinite( ) &&
                         ( enuToTerrestrial->transpose( ) * *enuToTerrestrial - Eigen::Matrix3d::Identity( ) ).norm( ) < 1e-12 &&
                         std::abs( enuToTerrestrial->determinant( ) - 1.0 ) < 1e-12,
                 "local basis must be a constant proper rotation" );
    std::vector< StationHistoryLinearSolution > linear;
    std::set< std::string > identifiers;
    std::vector< std::pair< std::string, StationHistoryInterval > > velocityIntervals;
    std::set< std::string > velocityIdentifiers;
    for( const auto& record : discontinuities.discontinuities )
        if( record.kind == "V" && sameSite( record.key.site, site ) )
        {
            require( velocityIdentifiers.insert( record.key.solutionId ).second, "repeated V interval identity" );
            velocityIntervals.emplace_back( record.key.solutionId,
                                            StationHistoryInterval{ intervalEnd( record.start, true, discontinuities ),
                                                                    intervalEnd( record.end, false, discontinuities ) } );
        }
    std::sort( velocityIntervals.begin( ), velocityIntervals.end( ), []( const auto& a, const auto& b ) {
        return !a.second.startTai ? b.second.startTai.has_value( ) : ( b.second.startTai && *a.second.startTai < *b.second.startTai );
    } );
    for( std::size_t i = 0; i < velocityIntervals.size( ); ++i )
    {
        const auto& interval = velocityIntervals[ i ].second;
        require( !interval.startTai || !interval.endTai || *interval.startTai < *interval.endTai, "empty or reversed V interval" );
        if( i )
            require( velocityIntervals[ i - 1 ].second.endTai && interval.startTai &&
                             *velocityIntervals[ i - 1 ].second.endTai == *interval.startTai,
                     "V intervals overlap or contain a gap" );
    }
    // A selector disambiguates SITE/ID metadata only. Parameter-to-monument
    // association must be warranted by the supplied companion products.
    auto selectedSolutions = solutions;
    if( domesSelector )
    {
        const auto range = selectedSolutions.sites.equal_range( site );
        for( auto it = range.first; it != range.second; )
            if( it->second.domesId != *domesSelector )
                it = selectedSolutions.sites.erase( it );
            else
                ++it;
    }
    for( const auto& key : getSinexCartesianSolutionKeys( solutions ) )
    {
        if( !sameSite( key.site, site ) ) continue;
        StationHistoryLinearSolution solution;
        solution.source = getSinexCartesianSolution( selectedSolutions, key, true );
        solution.referenceTime = sourceEpoch( solution.source.positionReferenceEpoch );
        for( int c = 0; c < 3; ++c )
        {
            solution.positionMetres[ c ] = solution.source.positionMetres[ c ];
            solution.velocityMetresPerJulianYear[ c ] = *solution.source.velocityMetresPerSecond[ c ] * secondsPerJulianYear;
        }
        std::size_t count = 0;
        for( const auto& record : discontinuities.discontinuities )
            if( record.kind == "P" && sameKey( record.key, key ) )
            {
                ++count;
                solution.validity = { intervalEnd( record.start, true, discontinuities ),
                                      intervalEnd( record.end, false, discontinuities ) };
            }
        require( count == 1, "each Cartesian solution requires exactly one matching P validity interval" );
        count = 0;
        for( const auto& entry : velocityIntervals )
        {
            const auto& velocityInterval = entry.second;
            const bool containsStart = !velocityInterval.startTai ||
                    ( solution.validity.startTai && *solution.validity.startTai >= *velocityInterval.startTai );
            const bool containsEnd =
                    !velocityInterval.endTai || ( solution.validity.endTai && *solution.validity.endTai <= *velocityInterval.endTai );
            if( containsStart && containsEnd )
            {
                ++count;
                solution.velocitySolutionId = entry.first;
                solution.velocityValidity = velocityInterval;
            }
        }
        require( count == 1, "each P interval requires exactly one independent containing V interval" );
        for( const auto& previous : linear )
            if( previous.velocitySolutionId == solution.velocitySolutionId )
                require( ( previous.velocityMetresPerJulianYear - solution.velocityMetresPerJulianYear ).cwiseAbs( ).maxCoeff( ) <=
                                 velocityConsistencyToleranceMetresPerJulianYear,
                         "attached SSC velocities contradict their shared V interval" );
        identifiers.insert( key.solutionId );
        linear.push_back( std::move( solution ) );
    }
    for( const auto& velocityId : velocityIdentifiers )
        require( std::any_of( linear.begin( ), linear.end( ), [ & ]( const auto& s ) { return s.velocitySolutionId == velocityId; } ),
                 "V interval has no associated P solution" );
    for( const auto& record : discontinuities.discontinuities )
        if( sameSite( record.key.site, site ) && record.kind == "P" )
            require( identifiers.count( record.key.solutionId ), "position interval has no matching Cartesian solution" );

    std::vector< StationPostSeismicTerm > psd;
    if( postSeismic )
    {
        metadata.signals.push_back( "empirical post-seismic deformation" );
        std::set< std::size_t > paired;
        for( const auto& entry : postSeismic->parameters )
        {
            const auto& amplitude = entry.second;
            if( !sameSite( amplitude.key.site, site ) ) continue;
            if( amplitude.type.size( ) == 6 && ( amplitude.type.substr( 0, 5 ) == "TLOG_" || amplitude.type.substr( 0, 5 ) == "TEXP_" ) )
                continue;
            require(
                    amplitude.type.size( ) == 6 && ( amplitude.type.substr( 0, 5 ) == "ALOG_" || amplitude.type.substr( 0, 5 ) == "AEXP_" ),
                    "unsupported PSD parameter" );
            const auto partner = postSeismic->parameters.find( amplitude.index + 1 );
            require( partner != postSeismic->parameters.end( ), "orphan PSD amplitude parameter" );
            const auto& relaxation = partner->second;
            require( relaxation.type == "T" + amplitude.type.substr( 1 ) && sameKey( relaxation.key, amplitude.key ) &&
                             relaxation.referenceEpoch.civilUtc && amplitude.referenceEpoch.civilUtc &&
                             *relaxation.referenceEpoch.civilUtc == *amplitude.referenceEpoch.civilUtc,
                     "PSD source-index partner has inconsistent family/component/identity/event" );
            require( relaxation.originalUnit == "m" || relaxation.originalUnit == "y", "unsupported IGN PSD relaxation unit" );
            StationPostSeismicTerm term;
            term.function = amplitude.type.substr( 1, 3 ) == "LOG" ? StationPostSeismicFunction::logarithmic
                                                                   : StationPostSeismicFunction::exponential;
            term.eventTime = sourceEpoch( amplitude.referenceEpoch );
            term.amplitudeMetres = metres( amplitude ) * componentVector( amplitude.type.back( ), enuToTerrestrial );
            term.relaxationJulianYears = relaxation.originalValue;
            term.amplitudeParameterIndex = amplitude.index;
            term.relaxationParameterIndex = relaxation.index;
            term.provenance = amplitude.provenance;
            psd.push_back( term );
            require( paired.insert( amplitude.index ).second && paired.insert( relaxation.index ).second, "PSD parameter used twice" );
        }
        for( const auto& entry : postSeismic->parameters )
            if( sameSite( entry.second.key.site, site ) ) require( paired.count( entry.first ), "orphan PSD relaxation parameter" );
    }
    std::vector< StationSeasonalTerm > seasons;
    if( seasonal )
    {
        require( !metadata.origin.empty( ), "seasonal product requires an explicit frame origin (for example CM or CF)" );
        metadata.signals.push_back( "empirical annual and semiannual station coordinates" );
        std::map< std::pair< std::string, int >, StationSeasonalTerm > grouped;
        std::map< std::pair< std::string, int >, std::set< std::string > > present;
        for( const auto& entry : seasonal->parameters )
        {
            const auto& p = entry.second;
            if( !sameSite( p.key.site, site ) ) continue;
            require( p.type.size( ) == 6 && p.type[ 0 ] == 'A' && ( p.type[ 1 ] == '1' || p.type[ 1 ] == '2' ) &&
                             ( p.type.substr( 2, 3 ) == "COS" || p.type.substr( 2, 3 ) == "SIN" ),
                     "unsupported seasonal parameter" );
            require( identifiers.count( p.key.solutionId ), "seasonal solution has no matching position validity" );
            const auto key = std::make_pair( p.key.solutionId, p.type[ 1 ] - '0' );
            require( present[ key ].insert( p.type ).second, "duplicate seasonal component" );
            auto& term = grouped[ key ];
            term.solutionId = key.first;
            term.cyclesPerCalendarYear = key.second;
            const Eigen::Vector3d component = metres( p ) * componentVector( p.type.back( ), enuToTerrestrial );
            if( p.type.substr( 2, 3 ) == "COS" )
                term.cosineMetres += component;
            else
                term.sineMetres += component;
            term.parameterIndices.push_back( p.index );
        }
        for( const auto& entry : grouped )
        {
            const auto& fields = present.at( entry.first );
            require( fields.size( ) == 6, "incomplete seasonal sine/cosine component set" );
            const char family = fields.begin( )->back( );
            const std::string components = std::string( "XYZ" ).find( family ) != std::string::npos ? "XYZ" : "ENU";
            for( const char c : components )
                for( const std::string function : { "COS", "SIN" } )
                    require( fields.count( "A" + std::to_string( entry.first.second ) + function + c ), "mixed seasonal component bases" );
            seasons.push_back( entry.second );
        }
        for( const auto& id : identifiers )
            require( grouped.count( { id, 1 } ) && grouped.count( { id, 2 } ),
                     "missing annual/semiannual coefficients for position solution" );
    }
    return std::make_shared< StationHistoryTrajectory >(
            std::move( linear ), std::move( psd ), std::move( seasons ), std::move( metadata ) );
}
}  // namespace ground_stations
}  // namespace tudat

/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/astro/ground_stations/stationCoordinateHistory.h"
#include "tudat/astro/reference_frames/terrestrialFrameTransformation.h"
#include "tudat/math/basic/mathematicalConstants.h"
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
constexpr double secondsPerJulianYear = 31557600.0;
void requireValidStationHistory( bool condition, const std::string& message )
{
    if( !condition )
    {
        throw std::runtime_error( "Station history: " + message );
    }
}
double elapsedJulianYears( const Time& tai, const Time& reference )
{
    return static_cast< double >( static_cast< long double >( tai - reference ) / secondsPerJulianYear );
}
std::string uniqueDomes( const SinexStationHistory& history,
                         const SinexSiteKey& site,
                         const std::optional< std::string >& selector = std::nullopt )
{
    const auto range = history.sites.equal_range( site );
    std::set< std::string > identifiers;
    for( auto it = range.first; it != range.second; ++it )
    {
        identifiers.insert( it->second.domesId );
    }
    if( selector )
    {
        requireValidStationHistory( !selector->empty( ) && identifiers.count( *selector ), history.source + ": selected DOMES is absent" );
        return *selector;
    }
    requireValidStationHistory( identifiers.size( ) == 1 && !identifiers.begin( )->empty( ),
                                history.source + ": station/point has missing or ambiguous DOMES association" );
    return *identifiers.begin( );
}
void checkMetadata( const SinexStationHistory& product, const StationHistoryMetadata& metadata )
{
    requireValidStationHistory(
            !product.frame.empty( ) && !product.release.empty( ) && !product.origin.empty( ) && !product.tideSystem.empty( ),
            product.source + ": explicit frame, release, origin and tide-system metadata required" );
    requireValidStationHistory( product.frame == metadata.frame && product.release == metadata.release &&
                                        product.origin == metadata.origin && product.tideSystem == metadata.tideSystem,
                                product.source + ": mixed frame, release, origin or tide-system metadata" );
}
std::optional< Time > intervalEnd( const SinexEpoch& epoch, bool start, const SinexStationHistory& history )
{
    if( epoch.meaning == ( start ? SinexEpochMeaning::openStart : SinexEpochMeaning::openEnd ) )
    {
        return std::nullopt;
    }
    const SinexEpoch* resolved = &epoch;
    if( epoch.meaning == ( start ? SinexEpochMeaning::headerStartDefault : SinexEpochMeaning::headerEndDefault ) )
    {
        resolved = start ? &history.headerStart : &history.headerEnd;
    }
    requireValidStationHistory( resolved->meaning == SinexEpochMeaning::specified && resolved->civilUtc.has_value( ),
                                history.source + ": unresolved station validity endpoint" );
    return convertSinexCivilUtcToTai( *resolved->civilUtc );
}
Time sourceEpoch( const SinexEpoch& epoch )
{
    requireValidStationHistory( epoch.meaning == SinexEpochMeaning::specified && epoch.civilUtc.has_value( ),
                                "unspecified model reference/event epoch" );
    return convertSinexCivilUtcToTai( *epoch.civilUtc );
}
Eigen::Vector3d componentVector( char component, const std::optional< Eigen::Matrix3d >& enuToTerrestrial )
{
    const auto xyz = std::string( "XYZ" ).find( component );
    if( xyz != std::string::npos )
    {
        return Eigen::Vector3d::Unit( xyz );
    }
    const auto enu = std::string( "ENU" ).find( component );
    requireValidStationHistory( enu != std::string::npos && enuToTerrestrial.has_value( ),
                                "ENU coefficients require an explicit constant local basis" );
    return enuToTerrestrial->col( enu );
}
double metres( const SinexParameterRecord& parameter )
{
    requireValidStationHistory( parameter.siUnit == "m" && parameter.valueSi.has_value( ),
                                "unsupported amplitude unit at parameter " + std::to_string( parameter.index ) );
    return *parameter.valueSi;
}
}  // namespace

//! Sort the linear solutions, which must have consecutive validity intervals.
StationCoordinateHistory::StationCoordinateHistory( std::vector< StationHistoryLinearSolution > solutions,
                                                    std::vector< StationPostSeismicTerm > postSeismicTerms,
                                                    std::vector< StationSeasonalTerm > seasonalTerms,
                                                    StationHistoryMetadata metadata ):
    solutions_( std::move( solutions ) ), postSeismicTerms_( std::move( postSeismicTerms ) ), seasonalTerms_( std::move( seasonalTerms ) ),
    metadata_( std::move( metadata ) )
{
    requireValidStationHistory( !solutions_.empty( ), "no linear position solutions supplied" );
    std::sort( solutions_.begin( ), solutions_.end( ), []( const auto& a, const auto& b ) {
        return !a.validity.startTai ? b.validity.startTai.has_value( )
                                    : ( b.validity.startTai && *a.validity.startTai < *b.validity.startTai );
    } );
    for( std::size_t i = 1; i < solutions_.size( ); ++i )
    {
        const auto& previous = solutions_[ i - 1 ].validity;
        const auto& interval = solutions_[ i ].validity;
        requireValidStationHistory( previous.endTai && interval.startTai && *previous.endTai == *interval.startTai,
                                    "position validity intervals overlap or contain a gap" );
    }
}

//! Evaluate owned linear, cumulative post-seismic and fitted seasonal records.
Eigen::Matrix< double, 6, 1 > StationCoordinateHistory::getState( const Time& tai ) const
{
    const StationHistoryLinearSolution* selected = nullptr;
    for( const auto& solution : solutions_ )
    {
        if( ( !solution.validity.startTai || tai >= *solution.validity.startTai ) &&
            ( !solution.validity.endTai || tai < *solution.validity.endTai ) )
        {
            selected = &solution;
        }
    }
    requireValidStationHistory( selected != nullptr, "requested epoch outside physical station validity" );
    const double yearsRate = 1.0 / secondsPerJulianYear;
    Eigen::Matrix< double, 6, 1 > result = Eigen::Matrix< double, 6, 1 >::Zero( );
    result.head< 3 >( ) =
            selected->positionMetres + selected->velocityMetresPerJulianYear * elapsedJulianYears( tai, selected->referenceTime );
    result.tail< 3 >( ) = selected->velocityMetresPerJulianYear * yearsRate;
    for( const auto& term : postSeismicTerms_ )
    {
        if( tai < term.eventTime )
        {
            continue;
        }
        const double dt = elapsedJulianYears( tai, term.eventTime );
        const double x = dt / term.relaxationJulianYears;
        if( term.function == StationPostSeismicFunction::logarithmic )
        {
            result.head< 3 >( ) += term.amplitudeMetres * std::log1p( x );
            result.tail< 3 >( ) += term.amplitudeMetres * ( yearsRate / ( term.relaxationJulianYears + dt ) );
        }
        else
        {
            result.head< 3 >( ) -= term.amplitudeMetres * std::expm1( -x );
            result.tail< 3 >( ) += term.amplitudeMetres * ( yearsRate / term.relaxationJulianYears * std::exp( -x ) );
        }
    }
    const long double decimalYear = getCalendarDecimalYear( tai );
    const int year = static_cast< int >( std::floor( decimalYear ) );
    const bool leapYear = year % 4 == 0 && ( year % 100 != 0 || year % 400 == 0 );
    const double decimalYearRate = 1.0 / ( ( leapYear ? 366.0 : 365.0 ) * 86400.0 );
    for( const auto& term : seasonalTerms_ )
    {
        if( term.solutionId != selected->source.key.solutionId )
        {
            continue;
        }
        const double omega = 2.0 * mathematical_constants::PI * term.cyclesPerCalendarYear;
        const double angle = omega * static_cast< double >( decimalYear - year );
        result.head< 3 >( ) += term.cosineMetres * std::cos( angle ) + term.sineMetres * std::sin( angle );
        result.tail< 3 >( ) +=
                ( -term.cosineMetres * std::sin( angle ) + term.sineMetres * std::cos( angle ) ) * ( omega * decimalYearRate );
    }
    return result;
}

namespace
{
//! Release of a realization in the catalogue: the release of its row, or else the name of the realization.
std::string catalogueRelease( const input_output::TerrestrialFrameTransformationData& catalogue, const std::string& realization )
{
    if( realization == catalogue.referenceFrame.realization )
    {
        return catalogue.referenceFrame.release;
    }
    const auto& identity = catalogue.getParameters( realization ).targetFrame;
    return identity.release.empty( ) ? identity.realization : identity.release;
}
}  // namespace

//! Transform the coefficients of a station coordinate history to another frame realization.
std::shared_ptr< StationCoordinateHistory > transformStationCoordinateHistory(
        const std::shared_ptr< StationCoordinateHistory >& history,
        const input_output::TerrestrialFrameTransformationData& catalogue,
        input_output::TerrestrialFrameIdentity target,
        const std::string& targetOrigin,
        const std::string& targetTideSystem,
        const bool acknowledgeAlignedReleases )
{
    requireValidStationHistory( history != nullptr, "a station coordinate history is required" );
    const auto& input = history->getMetadata( );

    // Origin and tide system are not converted, so they must agree with the input.
    if( target.realization.empty( ) || target.release.empty( ) || input.origin.empty( ) || input.tideSystem.empty( ) ||
        targetOrigin != input.origin || targetTideSystem != input.tideSystem )
    {
        throw std::runtime_error( "Station frame transformation requires explicit matching origin/tide conventions and target identity" );
    }

    // A release that differs from that of the catalogue must be declared aligned by the user.
    const bool sameRealization = input.frame == target.realization;
    if( sameRealization )
    {
        if( input.release != target.release && !acknowledgeAlignedReleases )
        {
            throw std::runtime_error( "Same-frame station release change requires an explicit alignment declaration" );
        }
    }
    else if( ( input.release != catalogueRelease( catalogue, input.frame ) ||
               target.release != catalogueRelease( catalogue, target.realization ) ) &&
             !acknowledgeAlignedReleases )
    {
        throw std::runtime_error( "Station releases differ from frame-catalogue identities; explicit alignment declaration required" );
    }

    std::vector< StationHistoryLinearSolution > solutions = history->getSolutions( );
    std::vector< StationPostSeismicTerm > postSeismicTerms = history->getPostSeismicTerms( );
    std::vector< StationSeasonalTerm > seasonalTerms = history->getSeasonalTerms( );
    if( !sameRealization )
    {
        const reference_frames::HelmertTransformationModel model( catalogue );

        // Each linear solution is transformed at its reference epoch; the parameter rates are per year (rate argument 1).
        for( auto& solution : solutions )
        {
            const auto transformation = model.getTransformation(
                    input.frame, target.realization, static_cast< double >( getCalendarDecimalYear( solution.referenceTime ) ), 1.0 );
            const Eigen::Vector3d position = solution.positionMetres;
            solution.positionMetres = transformation.translation + transformation.matrix * position;
            solution.velocityMetresPerJulianYear = transformation.translationRate + transformation.matrixRate * position +
                    transformation.matrix * solution.velocityMetresPerJulianYear;
            for( auto& term : seasonalTerms )
            {
                if( term.solutionId == solution.source.key.solutionId )
                {
                    term.cosineMetres = transformation.matrix * term.cosineMetres;
                    term.sineMetres = transformation.matrix * term.sineMetres;
                }
            }
        }

        // Post-seismic amplitudes are rotated and scaled with the matrix at the epoch of the earthquake.
        for( auto& term : postSeismicTerms )
        {
            term.amplitudeMetres =
                    model.getTransformation(
                                 input.frame, target.realization, static_cast< double >( getCalendarDecimalYear( term.eventTime ) ), 1.0 )
                            .matrix *
                    term.amplitudeMetres;
        }
    }

    StationHistoryMetadata metadata = input;
    metadata.frame = std::move( target.realization );
    metadata.release = std::move( target.release );
    return std::make_shared< StationCoordinateHistory >(
            std::move( solutions ), std::move( postSeismicTerms ), std::move( seasonalTerms ), std::move( metadata ) );
}

std::shared_ptr< StationCoordinateHistory > createIgnStationCoordinateHistory( const SinexStationHistory& solutions,
                                                                               const SinexStationHistory& discontinuities,
                                                                               const SinexSiteKey& site,
                                                                               const SinexStationHistory* postSeismic,
                                                                               const SinexStationHistory* seasonal,
                                                                               const std::optional< Eigen::Matrix3d >& enuToTerrestrial,
                                                                               const std::optional< std::string >& domesSelector )
{
    StationHistoryMetadata metadata;
    metadata.site = site;
    metadata.domesId = uniqueDomes( solutions, site, domesSelector );
    metadata.frame = solutions.frame;
    metadata.release = solutions.release;
    metadata.origin = solutions.origin;
    metadata.tideSystem = solutions.tideSystem;
    for( const auto* product : { &solutions, &discontinuities, postSeismic, seasonal } )
    {
        if( product )
        {
            checkMetadata( *product, metadata );
        }
    }

    // With a DOMES selector, only the selected monument is used to collect the solutions.
    auto selectedSolutions = solutions;
    if( domesSelector )
    {
        const auto range = selectedSolutions.sites.equal_range( site );
        for( auto it = range.first; it != range.second; )
        {
            it = it->second.domesId != *domesSelector ? selectedSolutions.sites.erase( it ) : std::next( it );
        }
    }

    // One linear solution per Cartesian solution of the station, valid in the interval of its position discontinuity record.
    std::vector< StationHistoryLinearSolution > linear;
    for( const auto& key : getSinexStationCoordinateSolutionKeys( solutions ) )
    {
        if( !( key.site == site ) )
        {
            continue;
        }
        StationHistoryLinearSolution solution;
        solution.source = getSinexStationCoordinateSolution( selectedSolutions, key, true );
        solution.referenceTime = sourceEpoch( solution.source.positionReferenceEpoch );
        for( int c = 0; c < 3; ++c )
        {
            solution.positionMetres[ c ] = solution.source.positionMetres[ c ];
            solution.velocityMetresPerJulianYear[ c ] = *solution.source.velocityMetresPerSecond[ c ] * secondsPerJulianYear;
        }
        const auto record = std::find_if( discontinuities.discontinuities.begin( ),
                                          discontinuities.discontinuities.end( ),
                                          [ & ]( const auto& entry ) { return entry.kind == "P" && entry.key == key; } );
        requireValidStationHistory( record != discontinuities.discontinuities.end( ),
                                    "solution " + key.solutionId + " has no position validity interval" );
        solution.validity = { intervalEnd( record->start, true, discontinuities ), intervalEnd( record->end, false, discontinuities ) };
        linear.push_back( std::move( solution ) );
    }

    // Post-seismic terms: each amplitude parameter (ALOG_/AEXP_) is followed by its relaxation time (TLOG_/TEXP_, in
    // Julian years, although IGN labels the unit m).
    std::vector< StationPostSeismicTerm > psd;
    for( const auto& entry : postSeismic ? postSeismic->parameters : std::map< std::size_t, SinexParameterRecord >( ) )
    {
        const auto& amplitude = entry.second;
        if( !( amplitude.key.site == site ) || amplitude.type.size( ) != 6 ||
            ( amplitude.type.substr( 0, 5 ) != "ALOG_" && amplitude.type.substr( 0, 5 ) != "AEXP_" ) )
        {
            continue;
        }
        const auto partner = postSeismic->parameters.find( amplitude.index + 1 );
        requireValidStationHistory( partner != postSeismic->parameters.end( ) && partner->second.type == "T" + amplitude.type.substr( 1 ),
                                    "post-seismic amplitude parameter " + std::to_string( amplitude.index ) + " has no relaxation time" );
        StationPostSeismicTerm term;
        term.function =
                amplitude.type.substr( 1, 3 ) == "LOG" ? StationPostSeismicFunction::logarithmic : StationPostSeismicFunction::exponential;
        term.eventTime = sourceEpoch( amplitude.referenceEpoch );
        term.amplitudeMetres = metres( amplitude ) * componentVector( amplitude.type.back( ), enuToTerrestrial );
        term.relaxationJulianYears = partner->second.originalValue;
        psd.push_back( term );
    }

    // Seasonal terms: annual and semi-annual cosine and sine amplitudes (A1COS_, A2SIN_, ...) per position solution.
    std::map< std::pair< std::string, int >, StationSeasonalTerm > grouped;
    for( const auto& entry : seasonal ? seasonal->parameters : std::map< std::size_t, SinexParameterRecord >( ) )
    {
        const auto& p = entry.second;
        if( !( p.key.site == site ) )
        {
            continue;
        }
        requireValidStationHistory( p.type.size( ) == 6 && p.type[ 0 ] == 'A' && ( p.type[ 1 ] == '1' || p.type[ 1 ] == '2' ) &&
                                            ( p.type.substr( 2, 3 ) == "COS" || p.type.substr( 2, 3 ) == "SIN" ),
                                    "unsupported seasonal parameter " + p.type );
        auto& term = grouped[ { p.key.solutionId, p.type[ 1 ] - '0' } ];
        term.solutionId = p.key.solutionId;
        term.cyclesPerCalendarYear = p.type[ 1 ] - '0';
        ( p.type.substr( 2, 3 ) == "COS" ? term.cosineMetres : term.sineMetres ) +=
                metres( p ) * componentVector( p.type.back( ), enuToTerrestrial );
    }
    std::vector< StationSeasonalTerm > seasons;
    for( const auto& entry : grouped )
    {
        seasons.push_back( entry.second );
    }
    return std::make_shared< StationCoordinateHistory >(
            std::move( linear ), std::move( psd ), std::move( seasons ), std::move( metadata ) );
}
}  // namespace ground_stations
}  // namespace tudat

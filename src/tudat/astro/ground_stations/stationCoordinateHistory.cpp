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

//! Own and validate the published solution and correction records once.
StationCoordinateHistory::StationCoordinateHistory( std::vector< StationHistoryLinearSolution > solutions,
                                                    std::vector< StationPostSeismicTerm > postSeismicTerms,
                                                    std::vector< StationSeasonalTerm > seasonalTerms,
                                                    StationHistoryMetadata metadata ):
    solutions_( std::move( solutions ) ), postSeismicTerms_( std::move( postSeismicTerms ) ), seasonalTerms_( std::move( seasonalTerms ) ),
    metadata_( std::move( metadata ) )
{
    requireValidStationHistory( !solutions_.empty( ), "no linear position solutions supplied" );
    requireValidStationHistory( !metadata_.frame.empty( ) && !metadata_.release.empty( ), "explicit frame and release metadata required" );
    std::sort( solutions_.begin( ), solutions_.end( ), []( const auto& a, const auto& b ) {
        return !a.validity.startTai ? b.validity.startTai.has_value( )
                                    : ( b.validity.startTai && *a.validity.startTai < *b.validity.startTai );
    } );
    std::set< std::string > identifiers;
    for( std::size_t i = 0; i < solutions_.size( ); ++i )
    {
        const auto& solution = solutions_[ i ];
        const auto& interval = solution.validity;
        requireValidStationHistory( solution.positionMetres.allFinite( ) && solution.velocityMetresPerJulianYear.allFinite( ),
                                    "nonfinite linear solution" );
        requireValidStationHistory( solution.source.key.site == metadata_.site && !solution.source.key.solutionId.empty( ) &&
                                            identifiers.insert( solution.source.key.solutionId ).second,
                                    "inconsistent or repeated position solution identity" );
        for( const auto& boundary : { interval.startTai, interval.endTai } )
        {
            requireValidStationHistory( !boundary || std::isfinite( static_cast< long double >( *boundary ) ),
                                        "nonfinite validity endpoint" );
        }
        requireValidStationHistory( !interval.startTai || !interval.endTai || *interval.startTai < *interval.endTai,
                                    "empty or reversed validity interval" );
        if( i > 0 )
        {
            const auto& previous = solutions_[ i - 1 ].validity;
            requireValidStationHistory( previous.endTai && interval.startTai && *previous.endTai == *interval.startTai,
                                        "position validity intervals overlap or contain a gap" );
        }
    }
    for( const auto& term : postSeismicTerms_ )
    {
        requireValidStationHistory(
                term.amplitudeMetres.allFinite( ) && std::isfinite( term.relaxationJulianYears ) && term.relaxationJulianYears > 0.0,
                "post-seismic relaxation must be positive and finite" );
        requireValidStationHistory(
                term.function == StationPostSeismicFunction::logarithmic || term.function == StationPostSeismicFunction::exponential,
                "unsupported post-seismic function" );
    }
    std::set< std::pair< std::string, int > > seasonalKeys;
    for( const auto& term : seasonalTerms_ )
    {
        requireValidStationHistory( identifiers.count( term.solutionId ) &&
                                            ( term.cyclesPerCalendarYear == 1 || term.cyclesPerCalendarYear == 2 ) &&
                                            seasonalKeys.emplace( term.solutionId, term.cyclesPerCalendarYear ).second,
                                    "unknown or repeated seasonal position-solution/frequency identity" );
        requireValidStationHistory( term.cosineMetres.allFinite( ) && term.sineMetres.allFinite( ), "nonfinite seasonal coefficient" );
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
    metadata.sources.push_back( catalogue.provenance.source );
    metadata.signals.push_back( sameRealization ? "declared same-frame identity" : "terrestrial frame transformation" );
    return std::make_shared< StationCoordinateHistory >(
            std::move( solutions ), std::move( postSeismicTerms ), std::move( seasonalTerms ), std::move( metadata ) );
}

std::shared_ptr< StationCoordinateHistory > createIgnStationCoordinateHistory( const SinexStationHistory& solutions,
                                                                               const SinexStationHistory& discontinuities,
                                                                               const SinexSiteKey& site,
                                                                               const SinexStationHistory* postSeismic,
                                                                               const SinexStationHistory* seasonal,
                                                                               const std::optional< Eigen::Matrix3d >& enuToTerrestrial,
                                                                               const std::optional< std::string >& domesSelector,
                                                                               double velocityConsistencyToleranceMetresPerJulianYear )
{
    requireValidStationHistory(
            std::isfinite( velocityConsistencyToleranceMetresPerJulianYear ) && velocityConsistencyToleranceMetresPerJulianYear >= 0.0,
            "invalid declared velocity consistency tolerance" );
    StationHistoryMetadata metadata;
    metadata.site = site;
    metadata.domesId = uniqueDomes( solutions, site, domesSelector );
    metadata.frame = solutions.frame;
    metadata.release = solutions.release;
    metadata.origin = solutions.origin;
    metadata.tideSystem = solutions.tideSystem;
    metadata.signals = { "piecewise linear station coordinates" };
    for( const auto* product : { &solutions, &discontinuities, postSeismic, seasonal } )
    {
        if( !product )
        {
            continue;
        }
        checkMetadata( *product, metadata );
        metadata.sources.push_back( product->source );
        if( product == postSeismic || product == seasonal )
        {
            requireValidStationHistory( uniqueDomes( *product, site, domesSelector ) == metadata.domesId,
                                        product->source + ": station monument association differs" );
        }
    }
    if( enuToTerrestrial )
    {
        requireValidStationHistory(
                enuToTerrestrial->allFinite( ) &&
                        ( enuToTerrestrial->transpose( ) * *enuToTerrestrial - Eigen::Matrix3d::Identity( ) ).norm( ) < 1e-12 &&
                        std::abs( enuToTerrestrial->determinant( ) - 1.0 ) < 1e-12,
                "local basis must be a constant proper rotation" );
    }
    std::vector< StationHistoryLinearSolution > linear;
    std::set< std::string > identifiers;
    std::vector< std::pair< std::string, StationHistoryInterval > > velocityIntervals;
    std::set< std::string > velocityIdentifiers;
    for( const auto& record : discontinuities.discontinuities )
    {
        if( record.kind == "V" && record.key.site == site )
        {
            requireValidStationHistory( velocityIdentifiers.insert( record.key.solutionId ).second, "repeated V interval identity" );
            velocityIntervals.emplace_back( record.key.solutionId,
                                            StationHistoryInterval{ intervalEnd( record.start, true, discontinuities ),
                                                                    intervalEnd( record.end, false, discontinuities ) } );
        }
    }
    std::sort( velocityIntervals.begin( ), velocityIntervals.end( ), []( const auto& a, const auto& b ) {
        return !a.second.startTai ? b.second.startTai.has_value( ) : ( b.second.startTai && *a.second.startTai < *b.second.startTai );
    } );
    for( std::size_t i = 0; i < velocityIntervals.size( ); ++i )
    {
        const auto& interval = velocityIntervals[ i ].second;
        requireValidStationHistory( !interval.startTai || !interval.endTai || *interval.startTai < *interval.endTai,
                                    "empty or reversed V interval" );
        if( i )
        {
            requireValidStationHistory( velocityIntervals[ i - 1 ].second.endTai && interval.startTai &&
                                                *velocityIntervals[ i - 1 ].second.endTai == *interval.startTai,
                                        "V intervals overlap or contain a gap" );
        }
    }

    // A selector disambiguates SITE/ID metadata only. Parameter-to-monument
    // association must be warranted by the supplied companion products.
    auto selectedSolutions = solutions;
    if( domesSelector )
    {
        const auto range = selectedSolutions.sites.equal_range( site );
        for( auto it = range.first; it != range.second; )
        {
            if( it->second.domesId != *domesSelector )
            {
                it = selectedSolutions.sites.erase( it );
            }
            else
            {
                ++it;
            }
        }
    }
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
        std::size_t count = 0;
        for( const auto& record : discontinuities.discontinuities )
        {
            if( record.kind == "P" && record.key == key )
            {
                ++count;
                solution.validity = { intervalEnd( record.start, true, discontinuities ),
                                      intervalEnd( record.end, false, discontinuities ) };
            }
        }
        requireValidStationHistory( count == 1, "each Cartesian solution requires exactly one matching P validity interval" );
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
        requireValidStationHistory( count == 1, "each P interval requires exactly one independent containing V interval" );
        for( const auto& previous : linear )
        {
            if( previous.velocitySolutionId == solution.velocitySolutionId )
            {
                requireValidStationHistory(
                        ( previous.velocityMetresPerJulianYear - solution.velocityMetresPerJulianYear ).cwiseAbs( ).maxCoeff( ) <=
                                velocityConsistencyToleranceMetresPerJulianYear,
                        "attached SSC velocities contradict their shared V interval" );
            }
        }
        identifiers.insert( key.solutionId );
        linear.push_back( std::move( solution ) );
    }
    for( const auto& velocityId : velocityIdentifiers )
    {
        requireValidStationHistory(
                std::any_of( linear.begin( ), linear.end( ), [ & ]( const auto& s ) { return s.velocitySolutionId == velocityId; } ),
                "V interval has no associated P solution" );
    }
    for( const auto& record : discontinuities.discontinuities )
    {
        if( record.key.site == site && record.kind == "P" )
        {
            requireValidStationHistory( identifiers.count( record.key.solutionId ),
                                        "position interval has no matching Cartesian solution" );
        }
    }

    std::vector< StationPostSeismicTerm > psd;
    if( postSeismic )
    {
        metadata.signals.push_back( "empirical post-seismic deformation" );
        std::set< std::size_t > paired;
        for( const auto& entry : postSeismic->parameters )
        {
            const auto& amplitude = entry.second;
            if( !( amplitude.key.site == site ) )
            {
                continue;
            }
            if( amplitude.type.size( ) == 6 && ( amplitude.type.substr( 0, 5 ) == "TLOG_" || amplitude.type.substr( 0, 5 ) == "TEXP_" ) )
            {
                continue;
            }
            requireValidStationHistory(
                    amplitude.type.size( ) == 6 && ( amplitude.type.substr( 0, 5 ) == "ALOG_" || amplitude.type.substr( 0, 5 ) == "AEXP_" ),
                    "unsupported PSD parameter" );
            const auto partner = postSeismic->parameters.find( amplitude.index + 1 );
            requireValidStationHistory( partner != postSeismic->parameters.end( ), "orphan PSD amplitude parameter" );
            const auto& relaxation = partner->second;
            requireValidStationHistory( relaxation.type == "T" + amplitude.type.substr( 1 ) && relaxation.key == amplitude.key &&
                                                relaxation.referenceEpoch.civilUtc && amplitude.referenceEpoch.civilUtc &&
                                                *relaxation.referenceEpoch.civilUtc == *amplitude.referenceEpoch.civilUtc,
                                        "PSD source-index partner has inconsistent family/component/identity/event" );
            requireValidStationHistory( relaxation.originalUnit == "m" || relaxation.originalUnit == "y",
                                        "unsupported IGN PSD relaxation unit" );
            StationPostSeismicTerm term;
            term.function = amplitude.type.substr( 1, 3 ) == "LOG" ? StationPostSeismicFunction::logarithmic
                                                                   : StationPostSeismicFunction::exponential;
            term.eventTime = sourceEpoch( amplitude.referenceEpoch );
            term.amplitudeMetres = metres( amplitude ) * componentVector( amplitude.type.back( ), enuToTerrestrial );
            term.relaxationJulianYears = relaxation.originalValue;
            term.amplitudeParameterIndex = amplitude.index;
            term.relaxationParameterIndex = relaxation.index;
            psd.push_back( term );
            requireValidStationHistory( paired.insert( amplitude.index ).second && paired.insert( relaxation.index ).second,
                                        "PSD parameter used twice" );
        }
        for( const auto& entry : postSeismic->parameters )
        {
            if( entry.second.key.site == site )
            {
                requireValidStationHistory( paired.count( entry.first ), "orphan PSD relaxation parameter" );
            }
        }
    }
    std::vector< StationSeasonalTerm > seasons;
    if( seasonal )
    {
        requireValidStationHistory( !metadata.origin.empty( ),
                                    "seasonal product requires an explicit frame origin (for example CM or CF)" );
        metadata.signals.push_back( "empirical annual and semiannual station coordinates" );
        std::map< std::pair< std::string, int >, StationSeasonalTerm > grouped;
        std::map< std::pair< std::string, int >, std::set< std::string > > present;
        for( const auto& entry : seasonal->parameters )
        {
            const auto& p = entry.second;
            if( !( p.key.site == site ) )
            {
                continue;
            }
            requireValidStationHistory( p.type.size( ) == 6 && p.type[ 0 ] == 'A' && ( p.type[ 1 ] == '1' || p.type[ 1 ] == '2' ) &&
                                                ( p.type.substr( 2, 3 ) == "COS" || p.type.substr( 2, 3 ) == "SIN" ),
                                        "unsupported seasonal parameter" );
            requireValidStationHistory( identifiers.count( p.key.solutionId ), "seasonal solution has no matching position validity" );
            const auto key = std::make_pair( p.key.solutionId, p.type[ 1 ] - '0' );
            requireValidStationHistory( present[ key ].insert( p.type ).second, "duplicate seasonal component" );
            auto& term = grouped[ key ];
            term.solutionId = key.first;
            term.cyclesPerCalendarYear = key.second;
            const Eigen::Vector3d component = metres( p ) * componentVector( p.type.back( ), enuToTerrestrial );
            if( p.type.substr( 2, 3 ) == "COS" )
            {
                term.cosineMetres += component;
            }
            else
            {
                term.sineMetres += component;
            }
            term.parameterIndices.push_back( p.index );
        }
        for( const auto& entry : grouped )
        {
            const auto& fields = present.at( entry.first );
            requireValidStationHistory( fields.size( ) == 6, "incomplete seasonal sine/cosine component set" );
            const char family = fields.begin( )->back( );
            const std::string components = std::string( "XYZ" ).find( family ) != std::string::npos ? "XYZ" : "ENU";
            for( const char c : components )
            {
                for( const std::string function : { "COS", "SIN" } )
                {
                    requireValidStationHistory( fields.count( "A" + std::to_string( entry.first.second ) + function + c ),
                                                "mixed seasonal component bases" );
                }
            }
            seasons.push_back( entry.second );
        }
        for( const auto& id : identifiers )
        {
            requireValidStationHistory( grouped.count( { id, 1 } ) && grouped.count( { id, 2 } ),
                                        "missing annual/semiannual coefficients for position solution" );
        }
    }
    return std::make_shared< StationCoordinateHistory >(
            std::move( linear ), std::move( psd ), std::move( seasons ), std::move( metadata ) );
}
}  // namespace ground_stations
}  // namespace tudat

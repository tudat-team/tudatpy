#include "tudat/io/readStationLoadingData.h"
#include "tudat/io/textParsing.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <locale>
#include <set>
#include <sstream>
#include <stdexcept>

namespace tudat
{
namespace input_output
{
namespace
{
const std::string signature = "EPHEDISP  Format version of 2005.06.30";
bool nextRecord( std::istream& stream, std::string& text )
{
    text.clear( );
    char character;
    while( stream.get( character ) )
    {
        if( character == '\n' || character == '\r' )
        {
            if( character == '\r' && stream.peek( ) == '\n' ) stream.get( );
            return true;
        }
        text.push_back( character );
    }
    return !text.empty( );
}
void failStationLoading( const std::string& source, const std::size_t line, const std::string& detail )
{
    throw std::runtime_error( "Station loading " + source + ( line ? " line " + std::to_string( line ) : "" ) + ": " + detail );
}
//! Parse a space-padded fixed-width field as an integer or (long) double of the requested type.
template< typename T >
T parseStationLoadingField( const std::string& field, const std::string& source, const std::size_t line )
{
    std::istringstream input( field );
    input.imbue( std::locale::classic( ) );
    T value{};
    std::string extra;
    if( !( input >> value ) || ( input >> extra ) ) failStationLoading( source, line, "invalid numeric field '" + field + "'" );
    if( !std::isfinite( static_cast< long double >( value ) ) ) failStationLoading( source, line, "nonfinite field" );
    return value;
}
void recordSize( const std::string& text, const std::size_t length, const std::string& source, const std::size_t line )
{
    if( text.size( ) < length || text.find_first_not_of( ' ', length ) != std::string::npos )
        failStationLoading( source, line, "truncated record or unexpected trailing field" );
}
void blanks( const std::string& text,
             const std::initializer_list< std::size_t > offsets,
             const std::string& source,
             const std::size_t line )
{
    for( const auto offset : offsets )
        if( text[ offset ] != ' ' ) failStationLoading( source, line, "invalid field delimiter" );
}
std::string label( const std::string& field, const std::string& source, const std::size_t line )
{
    const auto value = text_parsing::trim( field, " " );
    if( value.empty( ) || value.size( ) != field.find_last_not_of( ' ' ) + 1 || value.find( ' ' ) != std::string::npos )
        failStationLoading( source, line, "invalid source site identifier" );
    for( const unsigned char character : value )
        if( character < 33 ) failStationLoading( source, line, "invalid source site identifier character" );
    return value;
}
long double separation( const StationLoadingEpoch& first, const StationLoadingEpoch& second )
{
    return ( static_cast< long double >( second.modifiedJulianDay ) - first.modifiedJulianDay ) * 86400.0L + second.secondsOfDay -
            first.secondsOfDay;
}
StationLoadingEpoch indexedEpoch( const StationLoadingEpoch& begin, const long double offset, const std::string& source )
{
    const long double total = begin.secondsOfDay + offset;
    const long double days = std::floor( total / 86400.0L );
    const long double day = static_cast< long double >( begin.modifiedJulianDay ) + days;
    if( !std::isfinite( total ) || day < std::numeric_limits< int >::min( ) || day > std::numeric_limits< int >::max( ) )
        failStationLoading( source, 0, "epoch range overflow" );
    return { static_cast< int >( day ), total - days * 86400.0L };
}
}  // namespace

void validateStationLoadingData( const StationLoadingData& data )
{
    const auto& source = data.metadata.source;
    if( data.epochs.empty( ) || data.sites.empty( ) ) failStationLoading( source, 0, "empty epoch grid or sites" );
    if( data.timeScale != StationLoadingTimeScale::tai && data.timeScale != StationLoadingTimeScale::tt &&
        data.timeScale != StationLoadingTimeScale::utc && data.timeScale != StationLoadingTimeScale::tdb )
        failStationLoading( source, 0, "invalid declared source clock" );
    if( data.basis != StationLoadingBasis::geocentricRadialEastNorth && data.basis != StationLoadingBasis::cartesian &&
        data.basis != StationLoadingBasis::unspecified )
        failStationLoading( source, 0, "invalid declared source basis" );
    for( std::size_t i = 0; i < data.epochs.size( ); ++i )
    {
        const auto& epoch = data.epochs[ i ];
        if( !std::isfinite( epoch.secondsOfDay ) || epoch.secondsOfDay < 0.0L || epoch.secondsOfDay >= 86400.0L )
            failStationLoading( source, 0, "invalid split source epoch" );
        if( i && separation( data.epochs[ i - 1 ], epoch ) <= 0.0L ) failStationLoading( source, 0, "duplicate or unsorted source epochs" );
    }
    for( const auto declared : { data.declaredBeginEpoch, data.declaredEndEpoch } )
        if( declared &&
            ( !std::isfinite( declared->secondsOfDay ) || declared->secondsOfDay < 0.0L || declared->secondsOfDay >= 86400.0L ) )
            failStationLoading( source, 0, "invalid declared boundary epoch" );
    if( data.declaredBeginEpoch && separation( *data.declaredBeginEpoch, data.epochs.front( ) ) != 0.0L )
        failStationLoading( source, 0, "declared beginning disagrees with epoch grid" );
    const long double boundaryTolerance = .051L + data.epochs.size( ) * .5e-11L * 86400.0L;
    if( data.declaredEndEpoch && std::abs( separation( data.epochs.back( ), *data.declaredEndEpoch ) ) > boundaryTolerance )
        failStationLoading( source, 0, "declared end disagrees with epoch grid" );
    if( data.nominalSampleIntervalSeconds )
    {
        const auto interval = *data.nominalSampleIntervalSeconds;
        if( !std::isfinite( interval ) || interval < 0.0L || ( data.epochs.size( ) > 1 && interval == 0.0L ) )
            failStationLoading( source, 0, "invalid sampling interval" );
        for( std::size_t i = 1; i < data.epochs.size( ); ++i )
            if( std::abs( separation( data.epochs[ i - 1 ], data.epochs[ i ] ) - interval ) > 64.0L *
                        std::numeric_limits< long double >::epsilon( ) *
                        ( 1.0L + std::abs( data.epochs[ i - 1 ].secondsOfDay ) + std::abs( data.epochs[ i ].secondsOfDay ) +
                          std::abs( interval ) + std::abs( separation( data.epochs.front( ), data.epochs[ i ] ) ) ) )
                failStationLoading( source, 0, "epoch grid disagrees with declared sampling" );
    }
    std::set< StationLoadingComponent > components;
    if( data.metadata.referenceOrigin != StationLoadingOrigin::unspecified &&
        data.metadata.referenceOrigin != StationLoadingOrigin::wholeEarthCentreOfMass &&
        data.metadata.referenceOrigin != StationLoadingOrigin::solidEarthCentreOfMass &&
        data.metadata.referenceOrigin != StationLoadingOrigin::centreOfFigure &&
        data.metadata.referenceOrigin != StationLoadingOrigin::custom )
        failStationLoading( source, 0, "invalid reference origin" );
    for( const auto component : data.metadata.components )
    {
        if( component != StationLoadingComponent::atmosphere && component != StationLoadingComponent::hydrology &&
            component != StationLoadingComponent::ocean )
            failStationLoading( source, 0, "invalid physical component" );
        if( !components.insert( component ).second ) failStationLoading( source, 0, "duplicate physical component" );
    }
    for( const auto signal : { data.metadata.diurnal,
                               data.metadata.semidiurnal,
                               data.metadata.annual,
                               data.metadata.semiannual,
                               data.metadata.mean,
                               data.metadata.trend } )
        if( signal != StationLoadingSignalPresence::unknown && signal != StationLoadingSignalPresence::absent &&
            signal != StationLoadingSignalPresence::present )
            failStationLoading( source, 0, "invalid signal inventory" );
    std::set< std::string > harmonics;
    for( const auto& harmonic : data.metadata.removedHarmonics )
        if( harmonic.empty( ) || !harmonics.insert( harmonic ).second )
            failStationLoading( source, 0, "empty or duplicate removed harmonic" );
    for( const auto& entry : data.sites )
    {
        const auto& site = entry.second;
        if( entry.first.empty( ) || entry.first != site.identifier ) failStationLoading( source, site.sourceLine, "inconsistent site key" );
        for( const auto value : site.cartesianPosition )
            if( !std::isfinite( value ) ) failStationLoading( source, site.sourceLine, "nonfinite site location" );
        if( !std::isfinite( site.validityRadiusMetres ) || site.validityRadiusMetres < 0.0 )
            failStationLoading( source, site.sourceLine, "invalid validity radius" );
        std::optional< std::size_t > previous;
        for( const auto& sample : site.samples )
        {
            if( sample.first >= data.epochs.size( ) )
                failStationLoading( source, sample.second.sourceLine, "sample outside declared epoch grid" );
            if( previous && sample.first != *previous + 1 )
                failStationLoading( source, sample.second.sourceLine, "forbidden internal site gap" );
            for( const auto value : sample.second.displacement )
                if( !std::isfinite( value ) ) failStationLoading( source, sample.second.sourceLine, "nonfinite displacement" );
            previous = sample.first;
        }
    }
}

StationLoadingData readEphedispStationLoadingData( std::istream& stream, const StationLoadingMetadata& metadata )
{
    StationLoadingData data;
    data.metadata = metadata;
    if( data.metadata.source.empty( ) ) data.metadata.source = "<stream>";
    const auto& source = data.metadata.source;
    data.timeScale = StationLoadingTimeScale::tai;
    data.basis = StationLoadingBasis::geocentricRadialEastNorth;
    enum class Stage { header, counts, timeAndArea, sites, samples, trailer };
    Stage stage = Stage::header;
    std::size_t siteCount = 0, epochCount = 0, sampleCount = 0, readSamples = 0, lineNumber = 0;
    std::optional< StationLoadingEpoch > begin, end;
    std::optional< long double > interval;
    std::optional< double > radius;
    bool areaBeforeTime = false;
    std::optional< std::size_t > lastIndex;
    std::string text;
    while( nextRecord( stream, text ) )
    {
        ++lineNumber;
        if( !text.empty( ) && text.back( ) == '\r' ) text.pop_back( );
        if( stage == Stage::trailer ) failStationLoading( source, lineNumber, "record after trailer" );
        if( stage != Stage::header && !text.empty( ) && text.front( ) == '#' )
        {
            data.metadata.sourceHeaders.push_back( text );
            continue;
        }
        if( text.empty( ) ) failStationLoading( source, lineNumber, "empty record" );
        if( stage == Stage::header )
        {
            if( text != signature ) failStationLoading( source, lineNumber, "missing EPHEDISP signature" );
            stage = Stage::counts;
            continue;
        }
        if( stage == Stage::counts )
        {
            recordSize( text, 40, source, lineNumber );
            if( text.substr( 0, 8 ) != "P T 3 S " || text[ 19 ] != 'E' || text[ 28 ] != 'D' )
                failStationLoading( source, lineNumber, "invalid count record" );
            blanks( text, { 18, 20, 27, 29 }, source, lineNumber );
            const long long s = parseStationLoadingField< long long >( text.substr( 8, 10 ), source, lineNumber );
            const long long e = parseStationLoadingField< long long >( text.substr( 21, 6 ), source, lineNumber );
            const long long d = parseStationLoadingField< long long >( text.substr( 30, 10 ), source, lineNumber );
            if( s <= 0 || e <= 0 || d < 0 ) failStationLoading( source, lineNumber, "invalid declared counts" );
            siteCount = static_cast< std::size_t >( s );
            epochCount = static_cast< std::size_t >( e );
            sampleCount = static_cast< std::size_t >( d );
            if( static_cast< long double >( d ) > static_cast< long double >( s ) * e )
                failStationLoading( source, lineNumber, "displacement count exceeds grid" );
            stage = Stage::timeAndArea;
            continue;
        }
        if( text == signature )
        {
            if( stage != Stage::sites && stage != Stage::samples ) failStationLoading( source, lineNumber, "premature trailer" );
            stage = Stage::trailer;
            continue;
        }
        if( text.front( ) == 'A' || text.front( ) == 'T' )
        {
            if( stage != Stage::timeAndArea ) failStationLoading( source, lineNumber, "out-of-order time/area section" );
            if( text.front( ) == 'A' )
            {
                recordSize( text, 16, source, lineNumber );
                blanks( text, { 1 }, source, lineNumber );
                if( radius || ( ( begin || end || interval ) && !( begin && end && interval ) ) )
                    failStationLoading( source, lineNumber, "duplicate/interleaved area record" );
                radius = parseStationLoadingField< double >( text.substr( 2, 14 ), source, lineNumber );
                areaBeforeTime = !( begin || end || interval );
                if( *radius < 0.0 ) failStationLoading( source, lineNumber, "negative validity radius" );
            }
            else
            {
                const auto kind = text.substr( 0, 8 );
                if( radius && !areaBeforeTime ) failStationLoading( source, lineNumber, "time record after completed time/area section" );
                if( kind == "T sample" )
                {
                    recordSize( text, 26, source, lineNumber );
                    blanks( text, { 8, 9 }, source, lineNumber );
                    if( interval ) failStationLoading( source, lineNumber, "duplicate sample interval" );
                    interval = parseStationLoadingField< long double >( text.substr( 10, 16 ), source, lineNumber ) * 86400.0L;
                }
                else if( kind == "T begin " || kind == "T end   " )
                {
                    recordSize( text, 44, source, lineNumber );
                    blanks( text, { 8, 9, 15, 23, 24 }, source, lineNumber );
                    const StationLoadingEpoch epoch{ parseStationLoadingField< int >( text.substr( 10, 5 ), source, lineNumber ),
                                                     parseStationLoadingField< long double >( text.substr( 16, 7 ), source, lineNumber ) };
                    auto& target = kind == "T begin " ? begin : end;
                    if( target ) failStationLoading( source, lineNumber, "duplicate begin/end epoch" );
                    target = epoch;
                    ( kind == "T begin " ? data.declaredBeginEpoch : data.declaredEndEpoch ) = epoch;
                    ( kind == "T begin " ? data.informationalBeginEpoch : data.informationalEndEpoch ) = text.substr( 25, 19 );
                }
                else
                    failStationLoading( source, lineNumber, "unknown time record" );
            }
            if( begin && end && interval && radius )
            {
                if( *interval < 0 || ( epochCount > 1 && *interval == 0 ) || begin->secondsOfDay < 0 || begin->secondsOfDay >= 86400 ||
                    end->secondsOfDay < 0 || end->secondsOfDay >= 86400 )
                    failStationLoading( source, lineNumber, "invalid sampling or epoch" );
                // Source epochs have 0.1s precision and sample days have 1e-11 day precision.
                const long double tolerance = 0.051L + epochCount * 0.5e-11L * 86400.0L;
                if( std::abs( separation( *begin, *end ) - ( epochCount - 1 ) * *interval ) > tolerance )
                    failStationLoading( source, lineNumber, "declared end/count disagree with sampling" );
                data.nominalSampleIntervalSeconds = *interval;
                for( std::size_t i = 0; i < epochCount; ++i ) data.epochs.push_back( indexedEpoch( *begin, i * *interval, source ) );
                stage = Stage::sites;
            }
            continue;
        }
        if( text.front( ) == 'S' )
        {
            if( stage != Stage::sites ) failStationLoading( source, lineNumber, "out-of-order site section" );
            recordSize( text, 80, source, lineNumber );
            blanks( text, { 1, 2, 11, 12, 26, 40, 54, 55, 64, 73 }, source, lineNumber );
            StationLoadingSite site;
            site.identifier = label( text.substr( 3, 8 ), source, lineNumber );
            for( int i = 0; i < 3; ++i )
                site.cartesianPosition[ i ] = parseStationLoadingField< double >( text.substr( 13 + i * 14, 13 ), source, lineNumber );
            site.validityRadiusMetres = *radius;
            site.sourceLine = lineNumber;
            site.sourceRecord = text;
            site.informationalCoordinates = text.substr( 56, 24 );
            if( !data.sites.emplace( site.identifier, site ).second ) failStationLoading( source, lineNumber, "duplicate site" );
            if( data.sites.size( ) > siteCount ) failStationLoading( source, lineNumber, "excess site count" );
            continue;
        }
        if( text.front( ) == 'D' )
        {
            if( ( stage != Stage::sites && stage != Stage::samples ) || data.sites.size( ) != siteCount )
                failStationLoading( source, lineNumber, "out-of-order displacement section or missing sites" );
            stage = Stage::samples;
            recordSize( text, 80, source, lineNumber );
            blanks( text, { 1, 7, 8, 14, 22, 23, 43, 44, 53, 62, 71 }, source, lineNumber );
            const int index = parseStationLoadingField< int >( text.substr( 2, 5 ), source, lineNumber );
            if( index <= 0 || static_cast< std::size_t >( index ) > epochCount )
                failStationLoading( source, lineNumber, "sample index outside grid" );
            const auto zeroIndex = static_cast< std::size_t >( index - 1 );
            if( lastIndex && zeroIndex < *lastIndex ) failStationLoading( source, lineNumber, "displacements not sorted by epoch index" );
            lastIndex = zeroIndex;
            const auto identifier = label( text.substr( 45, 8 ), source, lineNumber );
            auto site = data.sites.find( identifier );
            if( site == data.sites.end( ) ) failStationLoading( source, lineNumber, "undeclared site" );
            StationLoadingSample sample;
            for( int i = 0; i < 3; ++i )
                sample.displacement[ i ] = parseStationLoadingField< double >( text.substr( 54 + i * 9, 8 ), source, lineNumber );
            sample.sourceLine = lineNumber;
            sample.sourceRecord = text;
            sample.informationalEpoch = text.substr( 9, 34 );
            if( !site->second.samples.emplace( zeroIndex, sample ).second )
                failStationLoading( source, lineNumber, "duplicate site/epoch sample" );
            if( ++readSamples > sampleCount ) failStationLoading( source, lineNumber, "excess displacement count" );
            continue;
        }
        failStationLoading( source, lineNumber, "unknown or out-of-order record" );
    }
    if( stream.bad( ) ) failStationLoading( source, lineNumber, "stream read failure" );
    if( stage != Stage::trailer ) failStationLoading( source, lineNumber, "missing trailer or incomplete sections" );
    if( data.sites.size( ) != siteCount || readSamples != sampleCount )
        failStationLoading( source, lineNumber, "declared record counts do not match" );
    validateStationLoadingData( data );
    return data;
}

StationLoadingData readEphedispStationLoadingData( const std::string& path, const StationLoadingMetadata& metadata )
{
    std::ifstream stream( path );
    if( !stream ) failStationLoading( path, 0, "cannot open file" );
    auto annotated = metadata;
    if( annotated.source.empty( ) ) annotated.source = path;
    return readEphedispStationLoadingData( stream, annotated );
}

StationLoadingMetadata imlsMerra2AtmosphericLoadingMetadata( const std::string& release )
{
    if( release.empty( ) ) throw std::invalid_argument( "IMLS MERRA2 loading profile requires a release identifier" );
    StationLoadingMetadata metadata;
    metadata.product = "IMLS atmospheric MERRA2 displacement series";
    metadata.release = release;
    metadata.sourceProfile = "IMLS MERRA2 atmosphere: pinned 2025-03-08 EPHEDISP description";
    metadata.referenceOrigin = StationLoadingOrigin::wholeEarthCentreOfMass;
    metadata.originDefinition = "IMLS whole-Earth centre of mass, including atmosphere";
    metadata.oceanResponse = "inverted barometer";
    metadata.components = { StationLoadingComponent::atmosphere };
    metadata.diurnal = metadata.semidiurnal = metadata.mean = StationLoadingSignalPresence::absent;
    metadata.annual = metadata.semiannual = metadata.trend = StationLoadingSignalPresence::present;
    metadata.removedHarmonics = { "PI1", "P1", "S1", "K1", "PSI1", "2T2", "T2", "S2" };
    return metadata;
}

}  // namespace input_output
}  // namespace tudat

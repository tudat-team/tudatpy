#include "tudat/io/readHighFrequencyEopData.h"

#include <cmath>
#include <fstream>
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
constexpr double microarcsecondsToRadians = 3.141592653589793238462643383279502884 / ( 180.0 * 3600.0 * 1.0e6 );

void fail( const std::string& source, const std::string& detail )
{
    throw std::runtime_error( "High-frequency EOP data " + source + ": " + detail );
}
}  // namespace

void validateHighFrequencyEopData( const HighFrequencyEopData& data )
{
    if( data.terms.empty( ) ) fail( data.source, "no tidal terms" );
    std::set< std::array< int, 6 > > multipliers;
    for( const auto& term : data.terms )
    {
        if( term.name.empty( ) || term.doodsonNumber.empty( ) ) fail( data.source, "missing tidal label" );
        if( !std::isfinite( term.periodDays ) || term.periodDays <= 0.0 ) fail( data.source, "invalid tidal period" );
        if( !multipliers.insert( term.argumentMultipliers ).second ) fail( data.source, "duplicate tidal argument multipliers" );
        for( unsigned int component = 0; component < 4; ++component )
        {
            if( !std::isfinite( term.sineAmplitudes[ component ] ) || !std::isfinite( term.cosineAmplitudes[ component ] ) )
                fail( data.source, "non-finite tidal amplitude" );
        }
    }
}

HighFrequencyEopData readHighFrequencyEopData( std::istream& stream, const std::string& source )
{
    HighFrequencyEopData data;
    data.source = source;
    std::string line;
    unsigned int lineNumber = 0;
    while( std::getline( stream, line ) )
    {
        ++lineNumber;
        const auto first = line.find_first_not_of( " \t\r" );
        if( first == std::string::npos ) continue;
        if( line[ first ] == '#' )
        {
            data.header.push_back( line );
            continue;
        }
        std::istringstream row( line );
        row.imbue( std::locale::classic( ) );
        HighFrequencyEopTerm term;
        row >> term.name;
        for( auto& multiplier : term.argumentMultipliers ) row >> multiplier;
        row >> term.doodsonNumber >> term.periodDays;
        for( unsigned int component = 0; component < 4; ++component )
        {
            row >> term.sineAmplitudes[ component ] >> term.cosineAmplitudes[ component ];
            const double factor = component < 2 ? microarcsecondsToRadians : 1.0e-6;
            term.sineAmplitudes[ component ] *= factor;
            term.cosineAmplitudes[ component ] *= factor;
        }
        std::string extra;
        if( row.fail( ) || ( row >> extra ) ) fail( source, "invalid 17-column record at line " + std::to_string( lineNumber ) );
        // Validate each line before adding it, to retain a useful location in errors.
        validateHighFrequencyEopData( { source + " line " + std::to_string( lineNumber ), {}, { term } } );
        data.terms.push_back( term );
    }
    if( stream.bad( ) ) fail( source, "stream read failed" );
    validateHighFrequencyEopData( data );
    return data;
}

HighFrequencyEopData readHighFrequencyEopData( const std::string& fileName )
{
    std::ifstream stream( fileName );
    if( !stream ) fail( fileName, "cannot open file" );
    return readHighFrequencyEopData( stream, fileName );
}

HighFrequencyEopData getDesaiSibois2016EopData( )
{
    // Factual coefficient table from the paper's ASCII supplement, mirrored by NASA IVS.
    // Exact source bytes SHA-256: 53398e40b70ccfa411d292e9a7f752a01cbac1be07d4e7fb18906d6481e57dd8.
    static const HighFrequencyEopData data = []( ) {
        const char* coefficients =
#include "desaiSibois2016EopData.inc"
                ;
        std::istringstream stream( coefficients );
        return readHighFrequencyEopData( stream,
                                         "https://ivscc.gsfc.nasa.gov/hfeop_wg/models/desai_model_jgrb51665-sup-0002-ds01.txt"
                                         "#sha256=53398e40b70ccfa411d292e9a7f752a01cbac1be07d4e7fb18906d6481e57dd8" );
    }( );
    return data;
}

}  // namespace input_output
}  // namespace tudat

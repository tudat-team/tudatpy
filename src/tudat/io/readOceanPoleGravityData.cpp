/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/io/readOceanPoleGravityData.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace tudat
{
namespace input_output
{
void validateOceanPoleGravityData( const OceanPoleGravityData& data )
{
    if( data.coefficients.empty( ) )
    {
        throw std::invalid_argument( "Empty ocean-pole gravity coefficient table" );
    }
    int maximumDegree = 0;
    for( const auto& item : data.coefficients )
    {
        const auto n = item.first.first, m = item.first.second;
        const auto& c = item.second;
        if( n < 1 || m < 0 || m > n || !std::isfinite( c.a.real( ) ) || !std::isfinite( c.a.imag( ) ) || !std::isfinite( c.b.real( ) ) ||
            !std::isfinite( c.b.imag( ) ) || ( m == 0 && c.b != std::complex< double >( 0.0 ) ) )
        {
            throw std::invalid_argument( "Invalid ocean-pole gravity coefficient at degree/order " + std::to_string( n ) + "/" +
                                         std::to_string( m ) );
        }
        maximumDegree = std::max( maximumDegree, n );
    }
    if( data.maximumDegree != maximumDegree )
    {
        throw std::invalid_argument( "Inconsistent ocean-pole gravity maximum degree" );
    }
}

OceanPoleGravityData readOceanPoleGravityData( std::istream& stream, const std::string& source )
{
    OceanPoleGravityData data;
    data.source = source;
    std::string line;
    std::size_t number = 0;
    while( std::getline( stream, line ) )
    {
        ++number;
        const auto start = line.find_first_not_of( " \t\r" );
        if( start == std::string::npos )
        {
            continue;
        }
        if( line[ start ] == '#' )
        {
            data.headerLines.push_back( line );
            continue;
        }
        if( data.coefficients.empty( ) )
        {
            std::istringstream header( line );
            std::string token, normalized;
            while( header >> token )
            {
                normalized += ( normalized.empty( ) ? "" : " " ) + token;
            }
            if( normalized == "n m Anm(Real) Bnm(Real) Anm(Imaginary) Bnm(Imaginary)" )
            {
                data.headerLines.push_back( line );
                continue;
            }
        }
        int n, m;
        double ar, br, ai, bi;
        std::string extra;
        std::istringstream row( line );
        if( !( row >> n >> m >> ar >> br >> ai >> bi ) || ( row >> extra ) )
        {
            throw std::invalid_argument( source + ":" + std::to_string( number ) + ": expected six ocean-pole gravity columns" );
        }
        if( !data.coefficients.emplace( std::make_pair( n, m ), OceanPoleGravityCoefficient{ { ar, ai }, { br, bi }, number } ).second )
        {
            throw std::invalid_argument( source + ":" + std::to_string( number ) + ": duplicate ocean-pole gravity degree/order" );
        }
        data.maximumDegree = std::max( data.maximumDegree, n );
    }
    if( stream.bad( ) )
    {
        throw std::runtime_error( "Cannot read ocean-pole gravity source " + source );
    }
    validateOceanPoleGravityData( data );
    return data;
}

OceanPoleGravityData readOceanPoleGravityData( const std::string& fileName )
{
    std::ifstream stream( fileName );
    if( !stream )
    {
        throw std::runtime_error( "Cannot open ocean-pole gravity source " + fileName );
    }
    return readOceanPoleGravityData( stream, fileName );
}
}  // namespace input_output
}  // namespace tudat

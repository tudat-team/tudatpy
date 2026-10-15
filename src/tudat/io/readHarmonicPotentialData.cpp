/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/io/readHarmonicPotentialData.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace tudat
{
namespace input_output
{
namespace
{
//! Published wave identity and inclusive degree/order extent for one identified provider profile.
struct ProfileWave {
    //! Canonical phase identifier and provider wave name.
    const char* doodson;

    const char* name;

    //! Inclusive extent of the published wave; zeros are permitted outside this extent.
    int minimumDegree, maximumDegree, maximumOrder;
};
#include "harmonicPotentialInventories.inc"

//! Return the pinned source inventory; generic IERS input has no inferred provider extent.
const std::vector< ProfileWave >& profileWaves( HarmonicPotentialProfile profile )
{
    static const std::vector< ProfileWave > unspecified;
    switch( profile )
    {
        case HarmonicPotentialProfile::iersPotential:
            return unspecified;
        case HarmonicPotentialProfile::fes2004:
            return fes2004Waves;
        case HarmonicPotentialProfile::fes2014bV1:
            return fes2014bV1Waves;
        case HarmonicPotentialProfile::costgEot11a:
            return costgEot11aWaves;
        default:
            throw std::invalid_argument( "Unknown harmonic potential profile" );
    }
}

//! Normalize an omitted leading zero and reject unverified Doodson character codes.
std::string canonicalDoodson( std::string identifier )
{
    if( identifier.size( ) == 6 )
    {
        identifier = "0" + identifier;
    }
    if( identifier.size( ) != 7 || identifier[ 3 ] != '.' )
    {
        throw std::invalid_argument( "Expected Doodson identifier ddd.ddd or dd.ddd" );
    }
    for( std::size_t i = 0; i < identifier.size( ); ++i )
    {
        if( i != 3 && ( identifier[ i ] < '0' || identifier[ i ] > '9' ) && identifier[ i ] != 'A' && identifier[ i ] != 'B' &&
            identifier[ i ] != 'Z' )
        {
            throw std::invalid_argument( "Invalid Doodson digit" );
        }
    }
    return identifier;
}
}  // namespace
//! Convert the provider identifier to the six phase multipliers without changing its phase convention.
std::array< int, 6 > parseDoodsonMultipliers( const std::string& identifier )
{
    const auto canonical = canonicalDoodson( identifier );

    std::array< int, 6 > multipliers{};

    // TiME22 revision 1.1 uses A=10, B=11 and Z=-1 before the usual five-digit offset.
    for( std::size_t i = 0; i < multipliers.size( ); ++i )
    {
        const char digit = canonical[ i < 3 ? i : i + 1 ];

        const int encoded = digit == 'A' ? 10 : digit == 'B' ? 11 : digit == 'Z' ? -1 : digit - '0';

        multipliers[ i ] = encoded - ( i == 0 ? 0 : 5 );
    }
    return multipliers;
}
//! Check sparse values and phase labels without inferring unpublished zeros or physical constants.
void validateHarmonicPotentialData( const HarmonicPotentialData& data )
{
    if( data.constituents.empty( ) || !std::isfinite( data.stokesUnitScale ) || data.stokesUnitScale <= 0.0 )
    {
        throw std::invalid_argument( "Missing harmonic potential inventory or Stokes units" );
    }
    std::set< std::string > ids;

    const auto& declared = profileWaves( data.profile );

    if( !declared.empty( ) && declared.size( ) != data.constituents.size( ) )
    {
        throw std::invalid_argument( "Harmonic potential inventory does not match the identified profile" );
    }
    if( !declared.empty( ) && data.stokesUnitScale != ( data.profile == HarmonicPotentialProfile::fes2004 ? 1e-12 : 1e-11 ) )
    {
        throw std::invalid_argument( "Harmonic potential units do not match the identified profile" );
    }
    for( std::size_t i = 0; i < data.constituents.size( ); ++i )
    {
        const auto& wave = data.constituents[ i ];

        if( canonicalDoodson( wave.doodson ) != wave.doodson || canonicalDoodson( wave.originalDoodson ) != wave.doodson ||
            parseDoodsonMultipliers( wave.doodson ) != wave.multipliers || wave.name.empty( ) || wave.coefficients.empty( ) ||
            !ids.insert( wave.doodson ).second )
        {
            throw std::invalid_argument( "Invalid or duplicate harmonic potential constituent" );
        }
        if( !declared.empty( ) &&
            ( wave.doodson != declared[ i ].doodson || wave.name != declared[ i ].name ||
              wave.minimumDegree != declared[ i ].minimumDegree || wave.maximumDegree != declared[ i ].maximumDegree ||
              wave.maximumOrder != declared[ i ].maximumOrder ) )
        {
            throw std::invalid_argument( "Harmonic potential ordered inventory/extent differs from the identified profile" );
        }
        for( const auto& item : wave.coefficients )
        {
            if( item.first.first < 0 || item.first.second < 0 || item.first.second > item.first.first )
            {
                throw std::invalid_argument( "Invalid harmonic potential degree/order" );
            }
            if( !declared.empty( ) &&
                ( item.first.first < wave.minimumDegree || item.first.first > wave.maximumDegree ||
                  item.first.second > wave.maximumOrder ) )
            {
                throw std::invalid_argument( "Harmonic potential coefficient lies outside its identified wave extent" );
            }
            for( double value : item.second.values )
            {
                if( !std::isfinite( value ) )
                {
                    throw std::invalid_argument( "Nonfinite harmonic potential coefficient" );
                }
            }
        }
    }
}
//! Read declared Stokes columns, keeping source order and checking the explicit provider profile.
HarmonicPotentialData readHarmonicPotentialData( std::istream& stream, const std::string& source, HarmonicPotentialProfile profile )
{
    HarmonicPotentialData data;

    data.source = source;
    data.profile = profile;
    std::map< std::string, std::size_t > indices;

    std::string line;

    std::size_t lineNumber = 0;

    bool columnsSeen = false, eotHeader = false, fes2004Header = false, fes2014Header = false;

    while( std::getline( stream, line ) )
    {
        ++lineNumber;
        std::istringstream tokens( line );

        std::string first;

        if( !( tokens >> first ) )
        {
            continue;
        }
        if( !columnsSeen )
        {
            data.headerLines.push_back( line );
            if( line.find( "Coefficients to compute variations in normalized Stokes coefficients (unit = 10^-11)" ) != std::string::npos )
            {
                if( data.stokesUnitScale != 0.0 )
                {
                    throw std::invalid_argument( "Duplicate harmonic potential units" );
                }
                data.stokesUnitScale = 1e-11;
            }
            if( line.find( "Coefficients to compute variations in normalized Stokes coefficients (unit = 10^-12)" ) != std::string::npos )
            {
                if( data.stokesUnitScale != 0.0 )
                {
                    throw std::invalid_argument( "Duplicate harmonic potential units" );
                }
                data.stokesUnitScale = 1e-12;
            }
            eotHeader |= line.find( "Ocean tide model: EOT11a up to (120,120) in cm" ) != std::string::npos;
            fes2004Header |= line.find( "Ocean tide model: FES2004 normalized model (fev. 2004) up to (100,100)" ) != std::string::npos;
            fes2014Header |= line.find( "Ocean tide model: FES2014b up to (100,100)" ) != std::string::npos;
            if( line.find( "in cm" ) != std::string::npos && ( profile != HarmonicPotentialProfile::costgEot11a || !eotHeader ) )
            {
                throw std::invalid_argument( "Water-height/contradictory units require the identified COST-G EOT11a potential profile" );
            }
            if( first == "Doodson" )
            {
                std::string normalized = first, token;

                while( tokens >> token )
                {
                    normalized += " " + token;
                }
                if( normalized != "Doodson Darw l m DelC+ DelS+ DelC- DelS-" || data.stokesUnitScale == 0.0 )
                {
                    throw std::invalid_argument( "Expected declared normalized Stokes units and IERS potential columns" );
                }
                columnsSeen = true;
            }
            continue;
        }
        if( first.front( ) == '#' )
        {
            data.headerLines.push_back( line );
            continue;
        }
        std::string name, extra;

        int n, m;

        std::array< double, 4 > values;

        if( !( tokens >> name >> n >> m >> values[ 0 ] >> values[ 1 ] >> values[ 2 ] >> values[ 3 ] ) || tokens >> extra )
        {
            throw std::invalid_argument( source + ":" + std::to_string( lineNumber ) + ": expected eight potential columns" );
        }
        // Normalize only the identifier spelling; phase and amplitude conversions belong to the consuming model.
        const auto id = canonicalDoodson( first );

        auto found = indices.find( id );

        if( found == indices.end( ) )
        {
            found = indices.emplace( id, data.constituents.size( ) ).first;
            data.constituents.push_back( { first, id, name, parseDoodsonMultipliers( id ), {} } );
        }
        auto& wave = data.constituents[ found->second ];

        if( wave.name != name ||
            !wave.coefficients.emplace( std::make_pair( n, m ), HarmonicPotentialCoefficient{ values, lineNumber } ).second )
        {
            throw std::invalid_argument( source + ":" + std::to_string( lineNumber ) +
                                         ": conflicting name or duplicate harmonic potential coefficient" );
        }
    }
    if( stream.bad( ) )
    {
        throw std::runtime_error( "Cannot read harmonic potential source " + source );
    }
    if( !columnsSeen || ( profile == HarmonicPotentialProfile::costgEot11a && ( !eotHeader || data.stokesUnitScale != 1e-11 ) ) ||
        ( profile == HarmonicPotentialProfile::fes2004 && ( !fes2004Header || data.stokesUnitScale != 1e-12 ) ) ||
        ( profile == HarmonicPotentialProfile::fes2014bV1 && ( !fes2014Header || data.stokesUnitScale != 1e-11 ) ) )
    {
        throw std::invalid_argument( "Harmonic potential header does not match the declared profile" );
    }
    const auto& declared = profileWaves( profile );

    if( !declared.empty( ) )
    {
        if( declared.size( ) != data.constituents.size( ) )
        {
            throw std::invalid_argument( "Harmonic potential inventory does not match the identified profile" );
        }
        for( std::size_t i = 0; i < declared.size( ); ++i )
        {
            data.constituents[ i ].minimumDegree = declared[ i ].minimumDegree;
            data.constituents[ i ].maximumDegree = declared[ i ].maximumDegree;
            data.constituents[ i ].maximumOrder = declared[ i ].maximumOrder;
        }
    }
    validateHarmonicPotentialData( data );
    return data;
}
//! Open the named potential file and preserve its path for source diagnostics.
HarmonicPotentialData readHarmonicPotentialData( const std::string& fileName, HarmonicPotentialProfile profile )
{
    std::ifstream stream( fileName );

    if( !stream )
    {
        throw std::runtime_error( "Cannot open harmonic potential source " + fileName );
    }
    return readHarmonicPotentialData( stream, fileName, profile );
}
}  // namespace input_output
}  // namespace tudat

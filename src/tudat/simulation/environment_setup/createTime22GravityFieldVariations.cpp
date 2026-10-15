/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/simulation/environment_setup/createTime22GravityFieldVariations.h"
#include "tudat/io/readIcgemGravityField.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace tudat
{
namespace simulation_setup
{
namespace
{
#include "time22Inventories.inc"

//! Read one of the profile's rectangular supplementary tables, rejecting missing or surplus values.
Eigen::MatrixXd readTime22Table( const std::filesystem::path& path, int rows, int columns )
{
    std::ifstream file( path );

    if( !file )
    {
        throw std::runtime_error( "Cannot open TiME22 table: " + path.string( ) );
    }
    Eigen::MatrixXd result( rows, columns );

    std::string line;

    for( int i = 0; i < rows; ++i )
    {
        if( !std::getline( file, line ) )
        {
            throw std::invalid_argument( "Missing TiME22 table row: " + path.string( ) );
        }
        std::istringstream row( line );

        for( int j = 0; j < columns; ++j )
        {
            if( !( row >> result( i, j ) ) || !std::isfinite( result( i, j ) ) )
            {
                throw std::invalid_argument( "Invalid TiME22 table value: " + path.string( ) );
            }
        }
        std::string extra;

        if( row >> extra )
        {
            throw std::invalid_argument( "Surplus TiME22 table columns: " + path.string( ) );
        }
    }
    while( std::getline( file, line ) )
    {
        if( line.find_first_not_of( " \t\r" ) != std::string::npos )
        {
            throw std::invalid_argument( "Surplus TiME22 table rows: " + path.string( ) );
        }
    }
    if( file.bad( ) )
    {
        throw std::runtime_error( "Cannot read TiME22 table: " + path.string( ) );
    }
    return result;
}

//! Check the product's physical reference and supplied full-potential, CF and phase declarations.
void validateTime22Field( const input_output::IcgemGravityFieldData& data,
                          const std::string& model,
                          const std::string& fileStem,
                          int maximumDegree )
{
    if( data.productType != "gravity_field" || data.modelName != model || data.gravitationalParameter != 3.986004418e14 ||
        data.referenceRadius != 6378136.6 || data.normalization != input_output::IcgemNormalization::fullyNormalized ||
        data.maximumDegree < maximumDegree )
    {
        throw std::invalid_argument( "Inconsistent TiME22 revision 1.1 field metadata: " + data.source );
    }
    bool loading = false, phase = false, centreOfFigure = false, identity = false;

    for( const auto& line : data.headerLines )
    {
        loading = loading || line.find( "Full tidal potential (includes loading potential)" ) != std::string::npos;
        phase = phase || line.find( "Doodson-Warburg phase correction is already applied." ) != std::string::npos;
        centreOfFigure = centreOfFigure || line.find( "Given in the Centre-of-figure (CF)" ) != std::string::npos;
        identity = identity || line.find( "--> " + fileStem ) != std::string::npos;
    }
    if( !loading || !phase || !centreOfFigure || !identity )
    {
        throw std::invalid_argument( "Missing TiME22 component, loading, CF or phase declaration: " + data.source );
    }
}
}  // namespace

//! Convert the declared paired product to final canonical coefficients and apply its ordered admittance once.
std::shared_ptr< GravityFieldVariationSettings > time22GravityFieldVariationSettings(
        const std::string& directory,
        gravitation::HarmonicGravityComponent component,
        int maximumDegree,
        int maximumOrder,
        bool applyAdmittance,
        const std::shared_ptr< ModelInterpolationSettings >& interpolation )
{
    using namespace gravitation;

    if( component != HarmonicGravityComponent::oceanMass && component != HarmonicGravityComponent::directAtmosphericMass )
    {
        throw std::invalid_argument( "TiME22 requires ocean mass or direct atmospheric mass" );
    }
    if( maximumDegree < 2 || maximumDegree > 360 || maximumOrder < 0 || maximumOrder > maximumDegree )
    {
        throw std::invalid_argument( "Invalid TiME22 Earth coefficient selection" );
    }
    const bool ocean = component == HarmonicGravityComponent::oceanMass;

    const std::string model = ocean ? "TiME22_OCN" : "TiME22_ATM";

    const auto& targets = ocean ? time22OceanTargets : time22AtmosphericTargets;

    const auto& inputs = ocean ? time22OceanInputs : time22AtmosphericTargets;

    const auto path = std::filesystem::path( directory ) / ( ocean ? "OCN" : "ATM" );

    // Match every supplied multiplier row to the pinned v1.1 inventory, retaining the supplied numerical phases.
    const auto table = readTime22Table( path / ( model + "_002doodson.txt" ), targets.size( ), 6 );

    std::map< std::string, std::array< int, 6 > > multipliers;

    for( std::size_t i = 0; i < targets.size( ); ++i )
    {
        const auto expected = input_output::parseDoodsonMultipliers( targets[ i ] );

        std::array< int, 6 > supplied{};

        for( int j = 0; j < 6; ++j )
        {
            if( table( i, j ) != expected[ j ] )
            {
                throw std::invalid_argument( "TiME22 revision 1.1 target multipliers or order do not match" );
            }
            supplied[ j ] = static_cast< int >( table( i, j ) );
        }
        multipliers.emplace( targets[ i ], supplied );
    }
    const auto publishedAdmittance = readTime22Table( path / ( model + "_003admittance.txt" ), inputs.size( ), targets.size( ) );

    if( !ocean && !publishedAdmittance.isIdentity( 0.0 ) )
    {
        throw std::invalid_argument( "TiME22 revision 1.1 atmospheric admittance must be the identity" );
    }
    std::ifstream fileList( path / ( model + "_001fileList.txt" ) );

    if( !fileList )
    {
        throw std::runtime_error( "Cannot open TiME22 ordered file list" );
    }
    std::vector< HarmonicGravityConstituent > constituents;

    std::string line;

    for( const auto& identifier : inputs )
    {
        std::string cosineFile, sineFile, extra;

        if( !std::getline( fileList, line ) )
        {
            throw std::invalid_argument( "Missing TiME22 input pair" );
        }
        std::istringstream row( line );

        if( !( row >> cosineFile >> sineFile ) || ( row >> extra ) || std::filesystem::path( cosineFile ).filename( ) != cosineFile ||
            cosineFile.compare( 0, model.size( ) + 9, model + "_" + identifier + "_" ) != 0 || cosineFile.size( ) < 8 ||
            cosineFile.substr( cosineFile.size( ) - 8 ) != "_cos.gfc" ||
            sineFile != cosineFile.substr( 0, cosineFile.size( ) - 7 ) + "sin.gfc" )
        {
            throw std::invalid_argument( "TiME22 input pair or order does not match revision 1.1" );
        }
        HarmonicGravityConstituent wave;

        wave.doodson = identifier;
        wave.multipliers = multipliers.at( identifier );
        wave.cCos = wave.cSin = wave.sCos = wave.sSin = Eigen::MatrixXd::Zero( maximumDegree + 1, maximumOrder + 1 );
        // Release each full raw field after copying the selected terms; no full-product raw cache is retained.
        for( bool cosine : { true, false } )
        {
            const auto& name = cosine ? cosineFile : sineFile;

            const auto data = input_output::readIcgemGravityField( ( path / name ).string( ) );

            validateTime22Field( data, model, std::filesystem::path( name ).stem( ).string( ), maximumDegree );
            for( int n = 2; n <= maximumDegree; ++n )
            {
                for( int m = 0; m <= std::min( n, maximumOrder ); ++m )
                {
                    const auto entry = data.coefficients.find( { n, m } );

                    if( entry == data.coefficients.end( ) )
                    {
                        throw std::invalid_argument( "Missing selected TiME22 coefficient: " + data.source );
                    }
                    ( cosine ? wave.cCos : wave.cSin )( n, m ) = entry->second.cosine;
                    ( cosine ? wave.sCos : wave.sSin )( n, m ) = entry->second.sine;
                }
            }
        }
        constituents.push_back( std::move( wave ) );
    }
    while( std::getline( fileList, line ) )
    {
        if( line.find_first_not_of( " \t\r" ) != std::string::npos )
        {
            throw std::invalid_argument( "Surplus TiME22 input pairs" );
        }
    }
    if( fileList.bad( ) )
    {
        throw std::runtime_error( "Cannot read TiME22 ordered file list" );
    }
    if( applyAdmittance )
    {
        // The published matrix is input-by-target; the existing coefficient map uses target-by-input.
        const HarmonicGravityAdmittance admittance{ inputs, targets, publishedAdmittance.transpose( ) };

        constituents = applyHarmonicGravityAdmittance( constituents, admittance );
        for( auto& wave : constituents )
        {
            wave.multipliers = multipliers.at( wave.doodson );
        }
    }
    const HarmonicGravityMetadata metadata{ component,
                                            3.986004418e14,
                                            6378136.6,
                                            true,
                                            "TiME22 revision 1.1 (2025-03-22), CF, loading and phase included",
                                            HarmonicGravityLoading::included };

    return std::make_shared< HarmonicGravityFieldVariationSettings >( constituents, metadata, maximumDegree, maximumOrder, interpolation );
}
}  // namespace simulation_setup
}  // namespace tudat

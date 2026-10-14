/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/astro/gravitation/harmonicGravityFieldVariations.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include <fstream>
#include <sstream>
using namespace tudat;
using namespace tudat::gravitation;
namespace
{
std::string fixture( const std::string& name, const std::string& directory = "harmonic_potential" )
{
    const std::string file = __FILE__;
    return file.substr( 0, file.find_last_of( "/" ) ) + "/../../../reference_files/" + directory + "/" + name;
}
HarmonicGravityMetadata metadata( )
{
    return { HarmonicGravityComponent::oceanMass, 1.0, 1.0, true, "controlled exact potential", HarmonicGravityLoading::included };
}
std::vector< std::vector< double > > rows( const std::string& path )
{
    std::ifstream stream( path );
    if( !stream )
    {
        throw std::runtime_error( "Missing harmonic gravity fixture" );
    }
    std::vector< std::vector< double > > result;
    std::string line;
    while( std::getline( stream, line ) )
    {
        std::istringstream input( line );
        std::vector< double > row;
        double value;
        while( input >> value )
        {
            row.push_back( value );
        }
        result.push_back( row );
    }
    return result;
}
}  // namespace
BOOST_AUTO_TEST_CASE( primaryPotentialSignsZonalConventionAndSingleGmRadiusConversion )
{
    std::istringstream stream(
            "Coefficients to compute variations in normalized Stokes coefficients (unit = 10^-12)\n"
            "Doodson Darw l m DelC+ DelS+ DelC- DelS-\n"
            "255.555 M2 2 0 5 6 0 0\n255.555 M2 2 1 7 8 9 10\n255.555 M2 2 2 1 2 3 4\n" );
    auto waves = convertHarmonicPotentialData( input_output::readHarmonicPotentialData( stream ), 2, 2 );
    Eigen::Vector6d arguments = Eigen::Vector6d::Zero( );
    auto callback = [ & ]( double ) { return arguments; };
    HarmonicGravityFieldVariations model( waves, metadata( ), callback, 1.0, 1.0, 2, 0, 2, 2 );
    auto a = model.calculateSphericalHarmonicsCorrections( 0 );
    BOOST_CHECK_SMALL( a.first( 0, 0 ) - 5e-12, 1e-25 );  // Published zonal amplitude is already doubled.
    BOOST_CHECK_SMALL( a.first( 0, 2 ) - 4e-12, 1e-25 );
    BOOST_CHECK_SMALL( a.second( 0, 2 ) + 2e-12, 1e-25 );
    BOOST_CHECK_EQUAL( a.second( 0, 0 ), 0.0 );
    arguments[ 0 ] = mathematical_constants::PI / 4;  // M2 theta=2*tau.
    a = model.calculateSphericalHarmonicsCorrections( 0 );
    BOOST_CHECK_SMALL( a.first( 0, 0 ) - 6e-12, 1e-25 );
    BOOST_CHECK_SMALL( a.first( 0, 2 ) - 6e-12, 1e-25 );
    BOOST_CHECK_SMALL( a.second( 0, 2 ) - 2e-12, 1e-25 );
    auto source = metadata( );
    source.gravitationalParameter = 2.0;
    source.referenceRadius = 2.0;
    HarmonicGravityFieldVariations scaled( waves, source, callback, 1.0, 1.0, 2, 1, 2, 2 );
    BOOST_CHECK_SMALL( scaled.calculateSphericalHarmonicsCorrections( 0 ).first( 0, 1 ) - 8 * a.first( 0, 2 ), 1e-25 );
    waves.front( ).cSin.setZero( );
    BOOST_CHECK_SMALL( ( model.calculateSphericalHarmonicsCorrections( 0 ).first - a.first ).norm( ), 1e-25 );
    source.gravitationalParameter = TUDAT_NAN;
    BOOST_CHECK_THROW( HarmonicGravityFieldVariations( waves, source, callback, 1.0, 1.0, 2, 0, 2, 2 ), std::invalid_argument );
}
BOOST_AUTO_TEST_CASE( namedWaveExtentsAllowOuterZerosAndRejectInternalHoles )
{
    auto data = input_output::readHarmonicPotentialData( fixture( "fes2014b-v1-degree3.potential" ),
                                                         input_output::HarmonicPotentialProfile::fes2014bV1 );
    const auto selected = convertHarmonicPotentialData( data, 3, 3, false, { "255.555", "055.575" } );
    BOOST_CHECK_EQUAL( selected[ 0 ].doodson, "255.555" );
    BOOST_CHECK_EQUAL( selected[ 1 ].doodson, "055.575" );
    BOOST_CHECK_SMALL( selected[ 1 ].cCos( 2, 0 ) + 0.06371e-11, 1e-25 );
    BOOST_CHECK_EQUAL( selected[ 1 ].cCos( 3, 0 ), 0.0 );
    for( auto& wave : data.constituents )
    {
        if( wave.doodson == "255.555" )
        {
            wave.coefficients.erase( { 3, 2 } );
        }
    }
    BOOST_CHECK_THROW( convertHarmonicPotentialData( data, 3, 3, true, { "255.555" } ), std::invalid_argument );
    BOOST_CHECK_THROW( convertHarmonicPotentialData( data, 3, 3, false, { "999.999" } ), std::invalid_argument );
    BOOST_CHECK_THROW( convertHarmonicPotentialData( data, 3, 3, false, { "055.575", "055.575" } ), std::invalid_argument );
}
BOOST_AUTO_TEST_CASE( admittanceUsesExactOrderedMajorsAndReplacesInventoryOnce )
{
    HarmonicGravityConstituent first;
    first.doodson = "255.555";
    first.multipliers = input_output::parseDoodsonMultipliers( first.doodson );
    first.cCos = first.cSin = first.sCos = first.sSin = Eigen::MatrixXd::Zero( 3, 3 );
    first.cCos( 2, 2 ) = 1e-12;
    auto second = first;
    second.doodson = "273.555";
    second.multipliers = input_output::parseDoodsonMultipliers( second.doodson );
    second.cCos( 2, 2 ) = 2e-12;
    HarmonicGravityAdmittance admittance{ { "255.555", "273.555" }, { "255.555", "273.555", "245.655" }, Eigen::MatrixXd( 3, 2 ) };
    admittance.matrix << 1, 0, 0, 1, 0.25, 0.75;
    const auto targets = applyHarmonicGravityAdmittance( { first, second }, admittance );
    BOOST_CHECK_EQUAL( targets.size( ), 3 );
    BOOST_CHECK_SMALL( targets[ 2 ].cCos( 2, 2 ) - 1.75e-12, 1e-25 );
    Eigen::Vector6d arguments = Eigen::Vector6d::Zero( );
    arguments[ 1 ] = mathematical_constants::PI;
    HarmonicGravityFieldVariations model( targets, metadata( ), [ & ]( double ) { return arguments; }, 1.0, 1.0, 2, 0, 2, 2 );
    // Major phases are zero and2*pi; N2 phase=-pi. Each target contributes once.
    BOOST_CHECK_SMALL( model.calculateSphericalHarmonicsCorrections( 0 ).first( 0, 2 ) - 1.25e-12, 1e-25 );
    std::swap( admittance.majorDoodsonIds[ 0 ], admittance.majorDoodsonIds[ 1 ] );
    BOOST_CHECK_THROW( applyHarmonicGravityAdmittance( { first, second }, admittance ), std::invalid_argument );
}
BOOST_AUTO_TEST_CASE( independentCostgEotM2Degree120ForcesWithSuppliedBenchmarkConstants )
{
    const auto data = input_output::readHarmonicPotentialData( fixture( "eot11a-degree3-with-full-m2.potential" ),
                                                               input_output::HarmonicPotentialProfile::costgEot11a );
    const auto waves = convertHarmonicPotentialData( data, 120, 120, false, { "255.555" } );
    const auto position = rows( fixture( "costg_00orbit_itrf.txt", "ocean_pole_gravity" ) ),
               rotation = rows( fixture( "costg_01earthRotation_rotaryMatrix.txt", "ocean_pole_gravity" ) ),
               arguments = rows( fixture( "costg_90doodsonArguments.txt" ) ),
               force = rows( fixture( "costg_10oceanTide_eot11a_M2_icrf.txt" ) );
    BOOST_REQUIRE_EQUAL( force.size( ), 4 );
    BOOST_REQUIRE_EQUAL( position.size( ), 4 );
    BOOST_REQUIRE_EQUAL( rotation.size( ), 4 );
    BOOST_REQUIRE_EQUAL( arguments.size( ), 4 );
    // Explicit matched benchmark convention, not inferred provider EOT11a reference constants.
    const double gm = 3.986004415e14, radius = 6378136.460;
    HarmonicGravityMetadata benchmark{
        HarmonicGravityComponent::oceanMass, gm, radius, true, "COST-G supplied EOT11a benchmark convention; EIGEN-6C4 GM/R",
        HarmonicGravityLoading::included
    };
    Eigen::Vector6d phase;
    HarmonicGravityFieldVariations model( waves, benchmark, [ & ]( double ) { return phase; }, gm, radius, 2, 0, 120, 120 );
    basic_mathematics::SphericalHarmonicsCache cache( 120, 120 );
    std::map< std::pair< int, int >, Eigen::Vector3d > perTerm;
    for( std::size_t i = 0; i < force.size( ); ++i )
    {
        for( int j = 0; j < 6; ++j )
        {
            phase[ j ] = arguments[ i ][ j + 1 ];
        }
        const auto block = model.calculateSphericalHarmonicsCorrections( 0 );
        Eigen::MatrixXd cosine = Eigen::MatrixXd::Zero( 121, 121 ), sine = cosine;
        cosine.block( 2, 0, 119, 121 ) = block.first;
        sine.block( 2, 0, 119, 121 ) = block.second;
        Eigen::Matrix3d matrix;
        for( int j = 0; j < 9; ++j )
        {
            matrix( j / 3, j % 3 ) = rotation[ i ][ j + 1 ];
        }
        const Eigen::Vector3d acceleration = computeGeodesyNormalizedGravitationalAccelerationSum(
                Eigen::Vector3d( position[ i ][ 1 ], position[ i ][ 2 ], position[ i ][ 3 ] ), gm, radius, cosine, sine, cache, perTerm );
        const Eigen::Vector3d residual =
                matrix.transpose( ) * acceleration - Eigen::Vector3d( force[ i ][ 1 ], force[ i ][ 2 ], force[ i ][ 3 ] );
        BOOST_TEST_MESSAGE( "COST-G M2 force residual " << residual.norm( ) << " m/s^2" );
        BOOST_CHECK_SMALL( residual.norm( ), 1e-11 );
    }
}

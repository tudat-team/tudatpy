#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/simulation/environment_setup/createGravityField.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/support/testFileUtilities.h"
#include <cstdio>
#include <fstream>
#include <sstream>

using namespace tudat;
namespace
{
input_output::IcgemGravityFieldData rawField( const std::string& norm )
{
    std::istringstream stream(
            "product_type gravity_field\nmodelname controlled\nearth_gravity_constant 3.986004415e14\n"
            "radius 6378136.46\nmax_degree 2\nerrors no\nnorm " +
            norm +
            "\ntide_system zero_tide\nend_of_head\n"
            "gfc 0 0 0.9 0\ngfc 1 0 0 0\ngfc 1 1 0 0\ngfc 2 0 -0.001 0\ngfc 2 1 0.000002 -0.000003\ngfc 2 2 0 0\n" );
    return input_output::readIcgemGravityField( stream );
}
}  // namespace

BOOST_AUTO_TEST_CASE( normalizationPresenceAndTideMetadataThroughStaticAndTimeDependentFactories )
{
    using basic_astrodynamics::TideSystem;
    const auto raw = rawField( "unnormalized" );
    const auto settings = simulation_setup::icgemSphericalHarmonicsGravitySettings( raw, "ITRF2020" );
    BOOST_CHECK( settings->getTideSystem( ) == TideSystem::zeroTide );
    BOOST_CHECK_EQUAL( settings->getCosineCoefficients( )( 0, 0 ), 0.9 );
    BOOST_CHECK_CLOSE_FRACTION( settings->getCosineCoefficients( )( 2, 0 ), -0.001 / std::sqrt( 5.0 ), 1e-15 );
    BOOST_CHECK_CLOSE_FRACTION( settings->getCosineCoefficients( )( 2, 1 ), 2e-6 * std::sqrt( 3.0 / 5.0 ), 1e-15 );
    BOOST_CHECK_CLOSE_FRACTION( settings->getSineCoefficients( )( 2, 1 ), -3e-6 * std::sqrt( 3.0 / 5.0 ), 1e-15 );
    BOOST_CHECK_EQUAL( raw.coefficients.at( { 2, 0 } ).cosine, -0.001 );
    const auto field = std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
            simulation_setup::createGravityFieldModel( settings, "Earth" ) );
    BOOST_REQUIRE( field );
    BOOST_CHECK( field->getTideSystem( ) == TideSystem::zeroTide );
    BOOST_CHECK_EQUAL( field->getFixedReferenceFrame( ), "ITRF2020" );
    BOOST_CHECK_EQUAL( field->getGravitationalParameter( ), raw.gravitationalParameter );
    BOOST_CHECK_EQUAL( field->getReferenceRadius( ), raw.referenceRadius );
    simulation_setup::SystemOfBodies bodies;
    bodies.createEmptyBody( "Earth" );
    settings->setCreateTimeDependentField( true );
    const auto varying = std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
            simulation_setup::createGravityFieldModel( settings, "Earth", bodies ) );
    BOOST_REQUIRE( varying );
    BOOST_CHECK( varying->getTideSystem( ) == TideSystem::zeroTide );
    const auto emptyVariations = std::make_shared< gravitation::GravityFieldVariationsSet >(
            std::vector< std::shared_ptr< gravitation::GravityFieldVariations > >{},
            std::vector< gravitation::BodyDeformationTypes >{},
            std::vector< std::string >{} );
    const gravitation::TimeDependentSphericalHarmonicsGravityField fullConstructor( raw.gravitationalParameter,
                                                                                    raw.referenceRadius,
                                                                                    settings->getCosineCoefficients( ),
                                                                                    settings->getSineCoefficients( ),
                                                                                    emptyVariations,
                                                                                    "ITRF2020",
                                                                                    TUDAT_NAN,
                                                                                    TideSystem::zeroTide );
    BOOST_CHECK( fullConstructor.getTideSystem( ) == TideSystem::zeroTide );
    BOOST_CHECK( gravitation::SphericalHarmonicsGravityField( 1.0, 1.0 ).getTideSystem( ) == TideSystem::unknown );
    auto sparse = raw;
    sparse.coefficients.erase( { 1, 0 } );
    sparse.coefficients.erase( { 1, 1 } );
    BOOST_CHECK_THROW( simulation_setup::icgemSphericalHarmonicsGravitySettings( sparse, "ITRF2020" ), std::runtime_error );
    const auto explicitZeros = simulation_setup::icgemSphericalHarmonicsGravitySettings( sparse, "ITRF2020", -1, -1, true );
    BOOST_CHECK_EQUAL( explicitZeros->getCosineCoefficients( )( 1, 0 ), 0.0 );
    BOOST_CHECK_EQUAL( sparse.coefficients.count( { 1, 0 } ), 0 );
    sparse.coefficients.erase( { 0, 0 } );
    BOOST_CHECK_THROW( simulation_setup::icgemSphericalHarmonicsGravitySettings( sparse, "ITRF2020", -1, -1, true ), std::runtime_error );
    const std::string source = __FILE__;
    const auto normalized = input_output::readIcgemGravityField( source.substr( 0, source.find_last_of( "/\\" ) ) +
                                                                 "/../../reference_files/icgem_gravity/gfz_rl061_ddk4_202104_degree2.gfc" );
    BOOST_CHECK_EQUAL(
            simulation_setup::icgemSphericalHarmonicsGravitySettings( normalized, "ITRF2020", 2, 2 )->getCosineCoefficients( )( 2, 0 ),
            -4.841653458641e-4 );

    // Even a static C00-only truncation must not silently discard a provider's temporal records.
    const auto temporal = input_output::readIcgemGravityField( source.substr( 0, source.find_last_of( "/\\" ) ) +
                                                               "/../../reference_files/icgem_gravity/eigen6s4_v2_c20_1950_2006.gfc" );

    BOOST_REQUIRE( !temporal.temporalCoefficients.empty( ) );
    BOOST_CHECK_THROW( simulation_setup::icgemSphericalHarmonicsGravitySettings( temporal, "ITRF2020", 0, 0 ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( numericReaderTruncatedOrderStillReadsLaterDegrees )
{
    const auto file = unit_tests::createTemporaryFilePath( "tudat-gravity", ".txt" );
    {
        std::ofstream stream( file );
        stream << "0 0 1 0\n1 0 0 0\n1 1 0 0\n2 0 -0.001 0\n2 1 2e-6 -3e-6\n2 2 99 99\n3 0 4e-6 0\n3 1 5e-6 6e-6\n";
    }
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > coefficients;
    simulation_setup::readGravityFieldFile( file, 3, 1, coefficients );
    std::remove( file.c_str( ) );
    BOOST_CHECK_EQUAL( coefficients.first.rows( ), 4 );
    BOOST_CHECK_EQUAL( coefficients.first.cols( ), 2 );
    BOOST_CHECK_EQUAL( coefficients.first( 2, 0 ), -0.001 );
    BOOST_CHECK_EQUAL( coefficients.first( 3, 0 ), 4e-6 );
    BOOST_CHECK_EQUAL( coefficients.first( 3, 1 ), 5e-6 );
    BOOST_CHECK_EQUAL( coefficients.second( 3, 1 ), 6e-6 );
}

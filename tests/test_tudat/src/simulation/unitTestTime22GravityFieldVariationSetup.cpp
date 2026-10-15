/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/simulation/environment_setup/createTime22GravityFieldVariations.h"
#include "tudat/io/readIcgemGravityField.h"
#include <cmath>
#include <fstream>
#include <sstream>

using namespace tudat;
using namespace tudat::simulation_setup;
using namespace tudat::gravitation;
namespace
{
//! Locate the small source extracts without requiring installed product resources.
std::string fixtureDirectory( )
{
    const std::string file = __FILE__;

    return file.substr( 0, file.find_last_of( "/" ) ) + "/../../reference_files/time22-v1_1";
}

//! Check all published phase rows against independently frozen provider contractions at four COST-G angle sets.
void checkProviderSynthesis( HarmonicGravityComponent component )
{
    const bool ocean = component == HarmonicGravityComponent::oceanMass;

    const std::string label = ocean ? "OCN" : "ATM";

    std::ifstream angleFile( fixtureDirectory( ) + "/../harmonic_potential/costg_90doodsonArguments.txt" );

    BOOST_REQUIRE( angleFile );
    std::vector< Eigen::Vector6d > arguments;

    std::string line;

    while( std::getline( angleFile, line ) )
    {
        std::istringstream row( line );

        double mjd;

        Eigen::Vector6d angle;

        BOOST_REQUIRE( row >> mjd );
        for( int i = 0; i < 6; ++i )
        {
            BOOST_REQUIRE( row >> angle[ i ] );
        }
        arguments.push_back( angle );
    }
    BOOST_REQUIRE_EQUAL( arguments.size( ), 4 );
    for( const bool admittance : { false, true } )
    {
        const auto base = time22GravityFieldVariationSettings( fixtureDirectory( ), component, 2, 2, admittance );

        const auto settings = std::dynamic_pointer_cast< HarmonicGravityFieldVariationSettings >( base );

        BOOST_REQUIRE( settings );
        BOOST_REQUIRE_EQUAL( settings->getConstituents( ).size( ), ocean ? ( admittance ? 435 : 57 ) : 16 );
        BOOST_CHECK_EQUAL( settings->getMetadata( ).gravitationalParameter, 3.986004418e14 );
        BOOST_CHECK_EQUAL( settings->getMetadata( ).referenceRadius, 6378136.6 );
        BOOST_CHECK( settings->getMetadata( ).loading == HarmonicGravityLoading::included );
        bool s1 = false, s5 = false, s6 = false;

        for( const auto& wave : settings->getConstituents( ) )
        {
            s1 = s1 || wave.doodson == ( ocean ? "164.555" : "164.556" );
            s5 = s5 || ( wave.doodson == "5A0.555" && wave.multipliers == std::array< int, 6 >{ 5, 5, -5, 0, 0, 0 } );
            s6 = s6 || ( wave.doodson == "6BZ.555" && wave.multipliers == std::array< int, 6 >{ 6, 6, -6, 0, 0, 0 } );
        }
        BOOST_CHECK( s1 && s5 && s6 );
        Eigen::Vector6d phaseArguments = Eigen::Vector6d::Zero( );

        HarmonicGravityFieldVariations model(
                settings->getConstituents( ),
                settings->getMetadata( ),
                [ &phaseArguments ]( double ) { return phaseArguments; },
                3.986004418e14,
                6378136.6,
                2,
                0,
                2,
                2 );

        std::ifstream expected( fixtureDirectory( ) + "/provider-reference.txt" );

        BOOST_REQUIRE( expected );
        std::getline( expected, line );
        int checked = 0;

        while( std::getline( expected, line ) )
        {
            std::istringstream row( line );

            std::string rowComponent;

            int useAdmittance, sample, degree, order;

            double cosine, sine, cosineBudget, sineBudget;

            BOOST_REQUIRE( row >> rowComponent >> useAdmittance >> sample >> degree >> order >> cosine >> sine >> cosineBudget >>
                           sineBudget );
            if( rowComponent != label || useAdmittance != admittance )
            {
                continue;
            }
            phaseArguments = arguments.at( sample );
            const auto actual = model.calculateSphericalHarmonicsCorrections( 0.0 );

            BOOST_CHECK_EQUAL( degree, 2 );
            BOOST_CHECK_LE( std::abs( actual.first( 0, order ) - cosine ), cosineBudget );
            BOOST_CHECK_LE( std::abs( actual.second( 0, order ) - sine ), sineBudget );
            ++checked;
        }
        BOOST_CHECK_EQUAL( checked, 12 );
    }
}
}  // namespace

BOOST_AUTO_TEST_CASE( atmosphericAllSixteenWavesIncludeExtendedCodesAndProviderPhase )
{
    checkProviderSynthesis( HarmonicGravityComponent::directAtmosphericMass );
    // Degree zero/one remain available in the neutral source record, while the Earth consuming block starts at two.
    const auto raw = input_output::readIcgemGravityField( fixtureDirectory( ) + "/ATM/TiME22_ATM_5A0.555_s5_cos.gfc" );

    BOOST_CHECK_NE( raw.coefficients.at( { 1, 1 } ).cosine, 0.0 );
    const auto settings = std::dynamic_pointer_cast< HarmonicGravityFieldVariationSettings >(
            time22GravityFieldVariationSettings( fixtureDirectory( ), HarmonicGravityComponent::directAtmosphericMass, 2, 2 ) );

    BOOST_CHECK_SMALL( settings->getConstituents( ).front( ).cCos.topRows( 2 ).norm( ), 1e-40 );
}

BOOST_AUTO_TEST_CASE( oceanAllFiftySevenInputsMapToFourHundredThirtyFiveTargetsOnce )
{
    checkProviderSynthesis( HarmonicGravityComponent::oceanMass );
}

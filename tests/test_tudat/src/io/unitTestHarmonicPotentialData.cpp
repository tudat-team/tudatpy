/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/io/readHarmonicPotentialData.h"
#include <sstream>
using namespace tudat::input_output;
namespace
{
std::string fixture( const std::string& name )
{
    const std::string file = __FILE__;
    return file.substr( 0, file.find_last_of( "/" ) ) + "/../../reference_files/harmonic_potential/" + name;
}
const std::string header =
        "Coefficients to compute variations in normalized Stokes coefficients (unit = 10^-12)\n"
        "Doodson Darw l m DelC+ DelS+ DelC- DelS-\n";
}  // namespace
BOOST_AUTO_TEST_CASE( publishedProfilesPreserveUnitsInventoryAndDeclaredWaveLimits )
{
    const auto fes2004 = readHarmonicPotentialData( fixture( "fes2004-degree3.potential" ), HarmonicPotentialProfile::fes2004 );
    const auto fes2014 = readHarmonicPotentialData( fixture( "fes2014b-v1-degree3.potential" ), HarmonicPotentialProfile::fes2014bV1 );
    const auto eot = readHarmonicPotentialData( fixture( "eot11a-degree3-with-full-m2.potential" ), HarmonicPotentialProfile::costgEot11a );
    BOOST_CHECK_EQUAL( fes2004.constituents.size( ), 18 );
    BOOST_CHECK_EQUAL( fes2014.constituents.size( ), 36 );
    BOOST_CHECK_EQUAL( eot.constituents.size( ), 18 );
    BOOST_CHECK_EQUAL( fes2004.stokesUnitScale, 1e-12 );
    BOOST_CHECK_EQUAL( fes2014.stokesUnitScale, 1e-11 );
    BOOST_CHECK_EQUAL( eot.stokesUnitScale, 1e-11 );
    BOOST_CHECK_EQUAL( fes2004.constituents.front( ).originalDoodson, "55.565" );
    BOOST_CHECK_EQUAL( fes2004.constituents.front( ).doodson, "055.565" );
    BOOST_CHECK_EQUAL( fes2004.constituents.front( ).maximumDegree, 2 );
    BOOST_CHECK_EQUAL( fes2014.constituents.front( ).maximumDegree, 100 );
    BOOST_CHECK_EQUAL( fes2014.constituents[ 1 ].coefficients.size( ), 1 );
    BOOST_CHECK_EQUAL( fes2014.constituents[ 1 ].maximumDegree, 2 );
    bool sa = false, s1 = false;
    for( const auto& wave : fes2014.constituents )
    {
        if( wave.name == "Sa" )
        {
            sa = wave.doodson == "056.555";
        }
        if( wave.name == "S1" )
        {
            s1 = wave.doodson == "164.555";
        }
    }
    BOOST_CHECK( sa );
    BOOST_CHECK( s1 );
    BOOST_CHECK_THROW( readHarmonicPotentialData( fixture( "eot11a-degree3-with-full-m2.potential" ) ), std::invalid_argument );
}
BOOST_AUTO_TEST_CASE( sixCharacterProviderIdsHaveCorrectDigitOffsets )
{
    const std::array< int, 6 > sa = { 0, 0, 1, 0, 0, -1 }, om1 = { 0, 0, 0, 0, 1, 0 };
    BOOST_CHECK( parseDoodsonMultipliers( "56.554" ) == sa );
    BOOST_CHECK( parseDoodsonMultipliers( "056.554" ) == sa );
    BOOST_CHECK( parseDoodsonMultipliers( "55.565" ) == om1 );
    const std::array< int, 6 > s5 = { 5, 5, -5, 0, 0, 0 }, s6 = { 6, 6, -6, 0, 0, 0 };

    BOOST_CHECK( parseDoodsonMultipliers( "5A0.555" ) == s5 );
    BOOST_CHECK( parseDoodsonMultipliers( "6BZ.555" ) == s6 );
    for( const auto id : { "055565", "055,565", "x55.565", "0055.565" } )
    {
        BOOST_CHECK_THROW( parseDoodsonMultipliers( id ), std::invalid_argument );
    }
}
BOOST_AUTO_TEST_CASE( ambiguousInvalidAndWaterHeightInputsReject )
{
    for( const auto records : { "55.565 Om1 2 0 1 2 3 4\n055.565 Om1 2 0 5 6 7 8\n",
                                "255.555 M2 2 3 1 2 3 4\n",
                                "255.555 M2 2 1 nan 2 3 4\n",
                                "255.555 M2 2 1 1 2 3\n",
                                "255.555 M2 2 1 1 2 3 4 extra\n" } )
    {
        std::istringstream stream( header + records );
        BOOST_CHECK_THROW( readHarmonicPotentialData( stream ), std::invalid_argument );
    }
    std::istringstream height( "Water heights in cm\nDoodson Darw l m DelC+ DelS+ DelC- DelS-\n255.555 M2 2 1 1 2 3 4\n" );
    BOOST_CHECK_THROW( readHarmonicPotentialData( height ), std::invalid_argument );
}

#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <cstdio>
#include <fstream>
#include "tudat/io/readVariousPdsFiles.h"
#include "tudat/support/testFileUtilities.h"

namespace tudat
{
namespace unit_tests
{

//! Read all three GRAIL mass formats while preserving their UTC or TDB epoch convention.
BOOST_AUTO_TEST_CASE( testGrailMassFormatsRetainTheirTimeConventions )
{
    const auto level0 = createTemporaryFilePath( "tudat-grail", "-mass0.txt" );

    const auto level1a = createTemporaryFilePath( "tudat-grail", "-mass1a.txt" );

    const auto level1b = createTemporaryFilePath( "tudat-grail", "-mass1b.txt" );

    // The same mission-era UTC epoch is written in the level-zero calendar and level-one-A numeric formats.
    {
        std::ofstream output( level0 );

        output << "header\nheader\nheader\nheader\n";
        output << "1   09/10/2011   13:08:50   13:08:52   0   0   201.5   0   0   0   0   0   0\n";
    }
    {
        std::ofstream output( level1a );

        for( int i = 0; i < 20; ++i )
        {
            output << "header\n";
        }
        output << "368932131.75 0.25 A 0 00000000 00000001 201.5\n";
    }
    {
        std::ofstream output( level1b );

        for( int i = 0; i < 20; ++i )
        {
            output << "header\n";
        }
        output << "368932132 0.25 A 0 00000000 00000001 201.5\n";
    }
    const auto calendar = input_output::grailMassLevel0FileReader( level0 );

    const auto numeric = input_output::grailMassLevel1FileReader( level1a, "1a" );

    const auto dynamical = input_output::grailMassLevel1FileReader( level1b, "1b" );

    std::remove( level0.c_str( ) );
    std::remove( level1a.c_str( ) );
    std::remove( level1b.c_str( ) );

    // The frozen prior-converter result retains the original 0.2 microsecond numerical budget.
    BOOST_REQUIRE_EQUAL( calendar.size( ), 1 );
    BOOST_REQUIRE_EQUAL( numeric.size( ), 1 );
    BOOST_REQUIRE_EQUAL( dynamical.size( ), 1 );
    BOOST_CHECK_SMALL( calendar.begin( )->first - 368932198.18251079, 2.0e-7 );
    BOOST_CHECK_SMALL( numeric.begin( )->first - 368932198.18251079, 2.0e-7 );

    // Level one B already supplies TDB; each reader retains the measured mass unchanged.
    BOOST_CHECK_EQUAL( dynamical.begin( )->first, 368932132.25 );
    BOOST_CHECK_EQUAL( calendar.begin( )->second, 201.5 );
    BOOST_CHECK_EQUAL( numeric.begin( )->second, 201.5 );
    BOOST_CHECK_EQUAL( dynamical.begin( )->second, 201.5 );
}

}  // namespace unit_tests
}  // namespace tudat

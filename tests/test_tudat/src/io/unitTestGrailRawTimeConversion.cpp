#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>

#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>

#include "tudat/io/readVariousPdsFiles.h"

namespace tudat
{
namespace unit_tests
{
namespace
{
struct TemporaryGrailFiles {
    TemporaryGrailFiles( )
    {
        directory = std::filesystem::temp_directory_path( ) /
                ( "tudat-grail-raw-" + std::to_string( std::chrono::steady_clock::now( ).time_since_epoch( ).count( ) ) );
        std::filesystem::create_directory( directory );
    }

    ~TemporaryGrailFiles( )
    {
        std::error_code error;
        std::filesystem::remove_all( directory, error );
    }

    std::filesystem::path directory;
};

// Frozen with TerrestrialTimeScaleConverter().getCurrentTime<double>(utc_scale,
// tdb_scale, utc, Eigen::Vector3d::Zero()) before removing that model dependency.
// Mission-era dates span the June 2012 leap-second transition. The tolerance is
// below a microsecond and accommodates the resolution of a mission-era double.
const std::array< std::array< double, 2 >, 4 > priorConverterEpochs = { {
        { 0.0, 64.183900714294609 },
        { 368932132.0, 368932198.18251079 },
        { 394372799.0, 394372865.18412083 },
        { 394372800.0, 394372867.18412083 },
} };
}  // namespace

BOOST_AUTO_TEST_CASE( testGrailLevel1aEpochsMatchPriorConverter )
{
    TemporaryGrailFiles temporary;
    const auto path = temporary.directory / "mass1a.txt";
    {
        std::ofstream output( path );
        for( int i = 0; i < 20; ++i )
        {
            output << "header\n";
        }
        output << std::setprecision( 17 );
        for( std::size_t i = 0; i < priorConverterEpochs.size( ); ++i )
        {
            output << priorConverterEpochs[ i ][ 0 ] - 0.25 << " 0.25 A 0 00000000 00000001 " << 200.0 + i << '\n';
        }
    }
    const auto masses = input_output::grailMassLevel1FileReader( path.string( ), "1a" );
    BOOST_REQUIRE_EQUAL( masses.size( ), priorConverterEpochs.size( ) );
    auto mass = masses.begin( );
    for( std::size_t i = 0; i < priorConverterEpochs.size( ); ++i, ++mass )
    {
        BOOST_CHECK_SMALL( mass->first - priorConverterEpochs[ i ][ 1 ], 2.0e-7 );
        BOOST_CHECK_EQUAL( mass->second, 200.0 + i );
    }
}

BOOST_AUTO_TEST_CASE( testGrailLevel0EpochMatchesPriorConverter )
{
    TemporaryGrailFiles temporary;
    const auto path = temporary.directory / "mass0.txt";
    {
        std::ofstream output( path );
        for( int i = 0; i < 4; ++i )
        {
            output << "header\n";
        }
        output << "1   09/10/2011   13:08:50   13:08:52   0   0   201.5   0   0   0   0   0   0\n";
    }
    const auto masses = input_output::grailMassLevel0FileReader( path.string( ) );
    BOOST_REQUIRE_EQUAL( masses.size( ), 1 );
    BOOST_CHECK_SMALL( masses.begin( )->first - priorConverterEpochs[ 1 ][ 1 ], 2.0e-7 );
    BOOST_CHECK_EQUAL( masses.begin( )->second, 201.5 );
}

BOOST_AUTO_TEST_CASE( testGrailLevel1bRetainsTdbEpoch )
{
    TemporaryGrailFiles temporary;
    const auto path = temporary.directory / "mass1b.txt";
    {
        std::ofstream output( path );
        for( int i = 0; i < 20; ++i )
        {
            output << "header\n";
        }
        output << "368932132 0.25 A 0 00000000 00000001 201.5\n";
    }
    const auto masses = input_output::grailMassLevel1FileReader( path.string( ), "1b" );
    BOOST_REQUIRE_EQUAL( masses.size( ), 1 );
    BOOST_CHECK_EQUAL( masses.begin( )->first, 368932132.25 );
    BOOST_CHECK_EQUAL( masses.begin( )->second, 201.5 );
}

}  // namespace unit_tests
}  // namespace tudat

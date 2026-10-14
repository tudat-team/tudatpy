/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/io/readOceanPoleGravityData.h"
#include <sstream>
using namespace tudat::input_output;
BOOST_AUTO_TEST_CASE( publishedCoefficientOrderAndSparsePresence )
{
    const std::string path = std::string( __FILE__ ).substr( 0, std::string( __FILE__ ).find_last_of( "/" ) ) +
            "/../../reference_files/ocean_pole_gravity/desaiscopolecoef_degree180.txt";
    const auto data = readOceanPoleGravityData( path );
    BOOST_CHECK_EQUAL( data.maximumDegree, 180 );
    BOOST_CHECK_EQUAL( data.coefficients.size( ), 16470 );
    BOOST_CHECK_EQUAL( data.headerLines.size( ), 1 );
    const auto c = data.coefficients.at( { 2, 2 } );
    BOOST_CHECK_EQUAL( c.a.real( ), 1.9102047023374e-2 );
    BOOST_CHECK_EQUAL( c.a.imag( ), -1.5123770169928e-2 );
    BOOST_CHECK_EQUAL( c.b.real( ), 1.1158297399424e-2 );
    BOOST_CHECK_EQUAL( c.b.imag( ), -2.4857839911518e-4 );
    BOOST_CHECK_EQUAL( c.sourceLine, 6 );
    std::istringstream sparse( "# separate gravity integrals\n2 1 0 0 0 0\n" );
    const auto zero = readOceanPoleGravityData( sparse );
    BOOST_CHECK_EQUAL( zero.coefficients.size( ), 1 );
    BOOST_CHECK_EQUAL( zero.coefficients.count( { 2, 0 } ), 0 );
}
BOOST_AUTO_TEST_CASE( ambiguousInvalidAndDisplacementRecordsReject )
{
    for( const auto text : { "2 1 1 2 3 4\n2 1 5 6 7 8\n",
                             "2 3 1 2 3 4\n",
                             "0 0 0 0 0 0\n",
                             "2 0 1 2 3 4\n",
                             "2 1 nan 2 3 4\n",
                             "2 1 1 2 3\n",
                             "2 1 1 2 3 4 5 6\n",
                             "arbitrary unknown header\n" } )
    {
        std::istringstream stream( text );
        BOOST_CHECK_THROW( readOceanPoleGravityData( stream ), std::invalid_argument );
    }
}

/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include <vector>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/InputOutput/missileDatcomReader.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_missile_datcom_data )

BOOST_AUTO_TEST_CASE( testMissileDatcomData )
{
    using namespace tudat;

    // Load missile Datcom data.
    std::string fileLocation = input_output::getTudatRootPath( )
            + "InputOutput/UnitTests/testFileMissileDatcomReader.dat";
    input_output::MissileDatcomReader missileDatcomReader( fileLocation );

    std::vector< double > missileDatcomData = missileDatcomReader.getMissileDatcomData( );
    double summation = 0.0;

    for ( unsigned int i = 0; i < missileDatcomData.size( ); i++ )
    {
        summation = summation + missileDatcomData[ i ];
    }

    // The input file consist of two times three sections. The first entry of each section is 1,
    // the next entry is the previous + 1. The first section is 144 long, the next 220 and the last
    // 400. When summed this gives the following result:
    const double expectedSummationResult = 229900.0;

    BOOST_CHECK_CLOSE_FRACTION( summation, expectedSummationResult,
                                std::numeric_limits< double >::epsilon( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

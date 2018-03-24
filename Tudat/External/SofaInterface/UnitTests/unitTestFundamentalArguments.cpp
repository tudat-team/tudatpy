/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <Eigen/LU>

#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/External/SofaInterface/fundamentalArguments.h"
#include "Tudat/External/SofaInterface/earthOrientation.h"


namespace tudat
{
namespace unit_tests
{

using namespace tudat::sofa_interface;

BOOST_AUTO_TEST_SUITE( test_sofa_fundamental_arguments )

//! Test calculation of Delaunay/Doodson arguments used for calculation of e.g. Earth-tide effects
BOOST_AUTO_TEST_CASE( testSofaFundamentalArguments )
{
    // Test whether Doodson <-> Delaunay argument conversion matrices are each other's inverse.
    Eigen::Matrix< double, 5, 5 > conversionMatrixComparison =
            delaunayToDoodsonArguments - doodsonToDelaunayArguments.inverse( );
    for( unsigned int i = 0; i < 5; i++ )
    {
        for( unsigned int j = 0; j < 5; j++ )
        {
            BOOST_CHECK_EQUAL( conversionMatrixComparison( i, j ), 0.0 );
        }
    }

    // Use test case from FUNDARG.F code (provided with IERS 2010 conventions) to calculate Delaunay arguments
    double testModifiedJulianDay1 = 54465.0;
    double testSecondsSinceJ2000 = ( testModifiedJulianDay1 -
            ( - basic_astrodynamics::JULIAN_DAY_AT_0_MJD + basic_astrodynamics::JULIAN_DAY_ON_J2000 ) ) *
            physical_constants::JULIAN_DAY;
    Eigen::Matrix< double, 5, 1 > expectedFundamentalArgumentValues;
    expectedFundamentalArgumentValues << 2.291187512612069099, 6.212931111003726414, 3.658025792050572989,
            4.554139562402433228, -0.5167379217231804489;

    // Calculate Delaunay arguments.
    Eigen::Matrix< double, 5, 1 > fundamentalArgumentValues =
            calculateDelaunayFundamentalArguments( testSecondsSinceJ2000 );

    // Compare against IERS results.
    for( unsigned int i = 0; i < 5; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( expectedFundamentalArgumentValues( i ) - fundamentalArgumentValues( i ) ), 2.0E-13  );
    }

    // Calculate Delaunay arguments with GMST.
    Eigen::Matrix< double, 6, 1 > fundamentalArgumentValuesWithGmst =
            calculateDelaunayFundamentalArgumentsWithGmst( testSecondsSinceJ2000 );
    for( unsigned int i = 0; i < 5; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( fundamentalArgumentValuesWithGmst( i + 1 ) - fundamentalArgumentValues( i ) ), 1.0E-15  );
    }

    // Manually compute GMST
    double expectedGmst = calculateGreenwichMeanSiderealTime(
                testSecondsSinceJ2000,
                convertTTtoUTC( testSecondsSinceJ2000 ),
                basic_astrodynamics::JULIAN_DAY_ON_J2000, basic_astrodynamics::iau_2006 );

    BOOST_CHECK_SMALL( std::fabs( expectedGmst + mathematical_constants::PI - fundamentalArgumentValuesWithGmst( 0 ) ), 1.0E-15  );

    // Calculate Doodson arguments directly and from Delaunay arguments and compare.
    Eigen::Matrix< double, 6, 1 > doodsonArguments = calculateDoodsonFundamentalArguments( testSecondsSinceJ2000 );
    Eigen::Matrix< double, 5, 1 > reconstructedFundamentalArguments =
            doodsonToDelaunayArguments * doodsonArguments.segment( 1, 5 );
    for( unsigned int i = 0; i < 5; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( fundamentalArgumentValues( i ) - reconstructedFundamentalArguments( i ) ), 1.0E-15  );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat




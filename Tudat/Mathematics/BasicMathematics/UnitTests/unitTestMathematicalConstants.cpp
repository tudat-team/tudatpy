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

#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{
namespace unit_tests
{

// Define Boost test suite.
BOOST_AUTO_TEST_SUITE( test_mathematical_constants )

//! Check correct pi using PI = circumference / diameter:
BOOST_AUTO_TEST_CASE( test_PI )
{    
    const double radius = 5.0;

    // Circumference circle with radius 5, 32 digits precision, see
    // http://www.wolframalpha.com/input/?i=Circumference+of+a+circle+with+radius+5
    // http://www.wolframalpha.com/input/?i=N[10*PI,66]
    double circumference = 31.4159265358979323846264338327950288419716939937510582097494459230; 
    BOOST_CHECK_CLOSE(  mathematical_constants::PI,
                        circumference / ( 2.0 * radius ) ,
                        std::numeric_limits< double >::epsilon( ) );
}

//! Check correct E using E Wolfram alpha as reference
BOOST_AUTO_TEST_CASE( test_E )
{    
    // Numerical value from:
    // http://www.wolframalpha.com/input/?i=e+72+digits
    BOOST_CHECK_CLOSE( mathematical_constants::E,
        2.71828182845904523536028747135266249775724709369995957496696762772407663, 
        std::numeric_limits< double >::epsilon( ) );
}

//! Check correct GOLDEN_RATIO using GOLDEN_RATIO Wolfram alpha as reference
BOOST_AUTO_TEST_CASE( test_GOLDEN_RATIO )
{    
    // Numerical value from:
    // http://www.wolframalpha.com/input/?i=golden+ratio+72+digits
    BOOST_CHECK_CLOSE(  mathematical_constants::GOLDEN_RATIO,
        1.618033988749894848204586834365638117720309179805762862135448622705260463, 
        std::numeric_limits< double >::epsilon( ) );
}

//! Check correct NAN using boost Floating Point Classification (fpclassify)
BOOST_AUTO_TEST_CASE( test_NAN )
{    
    // Numerical value from:
    // http://www.wolframalpha.com/input/?i=golden+ratio+72+digits
    BOOST_CHECK( boost::math::isnan( TUDAT_NAN ) );
}

BOOST_AUTO_TEST_CASE( test_TemplatedValues )
{
    double one = 1.0;
    long double longOne = 1.00000000000000000000L;

    BOOST_CHECK_CLOSE(  mathematical_constants::getFloatingInteger< double >( 1 ),
                        one, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE(  mathematical_constants::getFloatingInteger< long double >( 1 ),
                        longOne, std::numeric_limits< long  double >::epsilon( ) );

    BOOST_CHECK_SMALL(  mathematical_constants::getFloatingInteger< double >( 0 ),
                        std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL(  mathematical_constants::getFloatingInteger< long double >( 0 ),
                        std::numeric_limits< long  double >::epsilon( ) );

    double two = 2.0;
    long double longTwo = 2.00000000000000000000L;

    BOOST_CHECK_CLOSE(  mathematical_constants::getFloatingInteger< double >( 2 ),
                        two, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE(  mathematical_constants::getFloatingInteger< long double >( 2 ),
                        longTwo, std::numeric_limits< long  double >::epsilon( ) );

    BOOST_CHECK_CLOSE(  mathematical_constants::getFloatingFraction< double >( 1, 2 ),
                        mathematical_constants::getFloatingInteger< double >( 1 ) /
                        mathematical_constants::getFloatingInteger< double >( 2 ),
                        std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE(  mathematical_constants::getFloatingFraction< long double >( 1, 2 ),
                        mathematical_constants::getFloatingInteger< long double >( 1 ) /
                        mathematical_constants::getFloatingInteger< long double >( 2 )
                        , std::numeric_limits< long double >::epsilon( ) );


    BOOST_CHECK_CLOSE(  mathematical_constants::getPi< double >( ),
                        mathematical_constants::PI, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE(  mathematical_constants::getPi< long double >( ),
                        mathematical_constants::LONG_PI, std::numeric_limits< long  double >::epsilon( ) );

}

// Close Boost test suite.
BOOST_AUTO_TEST_SUITE_END( ) // End test_mathematical_constants

} // namespace unit_tests
} // namespace tudat

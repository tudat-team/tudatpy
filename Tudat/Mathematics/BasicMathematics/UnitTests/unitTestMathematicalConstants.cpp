/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      100910    J. Melman         First creation of code.
 *      110111    J. Melman         Adapted to the offical Tudat standards.
 *      110124    J. Melman         Further adapted to the offical Tudat standards.
 *      110201    J. Melman         Made the tests for obliquity and astronomical unit more
 *                                  accurate.
 *      120127    D. Dirkx          Moved to Tudat core.
 *      120127    K. Kumar          Transferred unit tests over to Boost unit test framework.
 *      120128    K. Kumar          Changed BOOST_CHECK to BOOST_CHECK_CLOSE_FRACTION for unit test
 *                                  comparisons.
 *      121205    K. Kumar          Updated license in file header.
 *
 *    References
 *
 *    Notes
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

// Close Boost test suite.
BOOST_AUTO_TEST_SUITE_END( ) // End test_mathematical_constants

} // namespace unit_tests
} // namespace tudat

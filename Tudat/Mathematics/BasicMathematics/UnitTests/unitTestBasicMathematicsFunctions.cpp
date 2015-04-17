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
 *      110207    B. Romgens        File created.
 *      110215    K. Kumar          Minor modifications to layout, comments
 *                                  and variable-naming.
 *      110411    K. Kumar          Added unit test for
 *                                  convertCartesianToSpherical( ) function.
 *      110701    K. Kumar          Updated failing tests with relative errors.
 *      110708    K. Kumar          Added unit tests for computeSampleMean( )
 *                                  and computeSampleVariance( ) functions.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      111111    K. Kumar          Strange error with convertCylindricalToCartesian function;
 *                                  achieved precision of results is less than machine precision,
 *                                  fixed by using slightly larger precision tolerance.
 *      120217    D. Dirkx          Bootified unit tests.
 *      120217    K. Kumar          Updated computeModuloForSignedValues() to computeModulo().
 *
 *    References
 *
 *    Notes
 *
 */

#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "TudatCore/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

namespace tudat
{
namespace unit_tests
{

using basic_mathematics::computeModulo;

BOOST_AUTO_TEST_SUITE( test_basic_mathematics )

//! Test if tudat modulo function is working.
BOOST_AUTO_TEST_CASE( testComputeModulo )
{
    // Test modulo function.
    // Test 1: Test 2.0 mod 2.0.
    // Test 2: Test 3.0 mod 2.5.
    // Test 3: Test 3.0 mod -2.5.
    // Test 4: Test -3.0 mod -2.5.
    // Test 5: Test -3.0 mod 2.5.

    const double machinePrecision = std::numeric_limits< double >::epsilon( );

    {
        const double resultUsingModuloFunction = computeModulo( 2.0, 2.0 );
        BOOST_CHECK_CLOSE_FRACTION( resultUsingModuloFunction, 0.0 , machinePrecision );
    }

    {
        const double resultUsingModuloFunction = computeModulo( 3.0, 2.5 );
        BOOST_CHECK_CLOSE_FRACTION( resultUsingModuloFunction, 0.5 , machinePrecision );
    }

    {
        const double resultUsingModuloFunction = computeModulo( 3.0, -2.5 );
        BOOST_CHECK_CLOSE_FRACTION( resultUsingModuloFunction, -2.0 , machinePrecision );
    }

    {
        const double resultUsingModuloFunction = computeModulo( -3.0, -2.5 );
        BOOST_CHECK_CLOSE_FRACTION( resultUsingModuloFunction, -0.5 , machinePrecision );
    }

    {
        const double resultUsingModuloFunction = computeModulo( -3.0, 2.5 );
        BOOST_CHECK_CLOSE_FRACTION( resultUsingModuloFunction, 2.0 , machinePrecision );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

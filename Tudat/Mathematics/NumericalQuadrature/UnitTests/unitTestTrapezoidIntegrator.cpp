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
 *      120202    K. Kumar          Moved unit tests from unitTestBasicMathematicsFunctions.h/.cpp;
 *                                  rewrote unit tests using Boost unit test framework.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/NumericalQuadrature/trapezoidIntegrator.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{
namespace unit_tests
{

using tudat::mathematical_constants::PI;

std::vector< double > linspace( double start, double end , int numberOfSamples )
{
    double spacing = ( end - start ) / ( static_cast< double >( numberOfSamples - 1 ) );

    std::vector< double > vector(0);
    for( int i = 0 ; i < numberOfSamples ; i++ )
    {
        vector.push_back( start + static_cast< double >(i) * spacing );
    }

    return vector;
}

BOOST_AUTO_TEST_SUITE( test_trapezoid_integrator )

//! Test if sample mean is computed correctly.
BOOST_AUTO_TEST_CASE( testIntegralSineFunction )
{
    // Generate independent variables
    int numberOfSamples = 1E4;

    std::vector< double > bounds(2);
    bounds[0] = 0 ;
    bounds[1] = PI;
    std::vector< double > independentVariables = linspace( bounds[0] , bounds[1] , numberOfSamples );

    std::vector< double > dependentVariables(0);
    for( unsigned int i = 0 ; i < independentVariables.size() ; i++ )
    {
        dependentVariables.push_back( std::sin( independentVariables[i] ) );
    }

    tudat::numerical_quadrature::TrapezoidNumericalIntegrator< double > integrator(
                                                                            independentVariables , dependentVariables );

    double expectedIntegral = 2.0;
    double computedIntegral = integrator.integrate();

    // Check if computed sample mean matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedIntegral, expectedIntegral, 1E-8 );
}

//! Test if sample mean is computed correctly.
BOOST_AUTO_TEST_CASE( testIntegralSineFunction2 )
{
    // Generate independent variables
    int numberOfSamples = 1E6;

    std::vector< double > bounds(2);
    bounds[0] = 0 ;
    bounds[1] = PI * 3.0 ;
    std::vector< double > independentVariables = linspace( bounds[0] , bounds[1] , numberOfSamples );

    std::vector< double > dependentVariables(0);
    for( unsigned int i = 0 ; i < independentVariables.size() ; i++ )
    {
        dependentVariables.push_back( std::sin( independentVariables[i] ) );
    }

    tudat::numerical_quadrature::TrapezoidNumericalIntegrator< double > integrator(
                                                                            independentVariables , dependentVariables );

    double expectedIntegral = 2.0;
    double computedIntegral = integrator.integrate();

    // Check if computed sample mean matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedIntegral, expectedIntegral, 1E-11 );
}

//! Test if sample mean is computed correctly.
BOOST_AUTO_TEST_CASE( testIntegralExpFunction )
{
    // Generate independent variables
    int numberOfSamples = 1E4;

    std::vector< double > bounds(2);
    bounds[0] = 0 ;
    bounds[1] = 2.0;
    std::vector< double > independentVariables = linspace( bounds[0] , bounds[1] , numberOfSamples );

    std::vector< double > dependentVariables(0);
    for( unsigned int i = 0 ; i < independentVariables.size() ; i++ )
    {
        dependentVariables.push_back( std::exp( independentVariables[i] ) );
    }

    tudat::numerical_quadrature::TrapezoidNumericalIntegrator< double > integrator(
                                                                            independentVariables , dependentVariables );

    double expectedIntegral = std::exp(2.0) - std::exp(0.0);
    double computedIntegral = integrator.integrate();

    // Check if computed sample mean matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedIntegral, expectedIntegral, 1E-8 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

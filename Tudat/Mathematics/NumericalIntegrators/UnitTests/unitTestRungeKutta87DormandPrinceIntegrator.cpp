/*    Copyright (c) 2010-2012, Delft University of Technology
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
 *      120323    K. Kumar          Created based off of unit test for RKF78 integrator; test data
 *                                  taken from (The Mathworks, 2012).
 *      120404    K. Kumar          Updated Matlab unit test by adding discrete-event data file.
 *
 *    References
 *      Montenbruck, O., Gill, E. Satellite Orbits: Models, Methods, Applications, Springer, 2005.
 *      The Mathworks, Inc. RKF78, Symbolic Math Toolbox, 2012.
 *
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h"
#include "Tudat/Mathematics/NumericalIntegrators/UnitTests/matlabNumericalIntegratorTest.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_runge_kutta_87_dormand_and_prince_integrator )

using numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd;
using numerical_integrators::RungeKuttaCoefficients;

//! Test Runge-Kutta 87 Dormand-Prince integrator using benchmark data from (The Mathworks, 2012).
BOOST_AUTO_TEST_CASE( testRungeKutta87DormandAndPrinceIntegratorUsingMatlabData )
{
    // Read in benchmark data (generated using Symbolic Math Toolbox in Matlab
    // (The Mathworks, 2012)). This data is generated using the DOPRI87 numerical integrator.
    //
    // It is important to note that the DOPRI87 coefficients are different from the coefficients
    // included in Tudat, which come from (Montenbruck and Gill, 2005). Nevertheless, this test
    // shows that they both produce near-identical results. If any problems are experienced with
    // this integrator, it might be wise to make use of the coefficients given in
    // (The Mathworks, 2012).
    std::string pathToBenchmarkDatafile = tudat::input_output::getTudatRootPath( )
            + "/Mathematics/NumericalIntegrators/UnitTests"
            + "/matlabOutputRungeKutta87DormandPrince.txt";

    // Read in discrete event benchmark data (generated using Symbolic Math Toolbox in Matlab
    // (The Mathworks, 2012)). This data is generated using the DOPRI87 numerical integrator.
    std::string pathToDiscreteEventBenchmarkDatafile = tudat::input_output::getTudatRootPath( )
            + "/Mathematics/NumericalIntegrators/UnitTests"
            + "/matlabOutputDiscreteEventRungeKutta87DormandPrince.txt";

    // Run Matlab integrator tests.
    matlab_numerical_integrator_tests::runMatlabNumericalIntegratorTests(
                pathToBenchmarkDatafile, std::numeric_limits< double >::epsilon( ), 1.0e-13,
                RungeKuttaCoefficients::get( RungeKuttaCoefficients::rungeKutta87DormandPrince ),
                pathToDiscreteEventBenchmarkDatafile );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120329    K. Kumar          Created based off of unit test for DOPRI87 integrator.
 *      120404    K. Kumar          Updated Matlab unit test by adding discrete-event data file.
 *
 *    References
 *      Montenbruck, O., Gill, E. Satellite Orbits: Models, Methods, Applications, Springer, 2005.
 *      The Mathworks, Inc. DOPRI56, Symbolic Math Toolbox, 2012.
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

BOOST_AUTO_TEST_SUITE( test_runge_kutta_fehlberg_56_integrator )

using mathematics::numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd;
using mathematics::numerical_integrators::RungeKuttaCoefficients;

//! Test Runge-Kutta-Fehlberg 56 integrator using benchmark data from (The Mathworks, 2012).
BOOST_AUTO_TEST_CASE( testRungeKuttaFehlberg56IntegratorUsingMatlabData )
{
    // Read in benchmark data (generated using Symbolic Math Toolbox in Matlab
    // (The Mathworks, 2012)). This data is generated using the DOPRI56 numerical integrator.
    //
    // It is important to note that the DOPRI56 coefficients are different from the RKF56
    // coefficients tested here, which come from (Montenbruck and Gill, 2005). Nevertheless, this
    // test shows that they both produce near-identical results. If any problems are experienced
    // with this integrator, it might be wise to find benchmark data specifically for the RKF56
    // integrator.
    std::string pathToBenchmarkDatafile = tudat::input_output::getTudatRootPath( )
            + "/Mathematics/NumericalIntegrators/UnitTests"
            + "/matlabOutputRungeKutta56DormandPrince.txt";

    // Read in discrete event benchmark data (generated using Symbolic Math Toolbox in Matlab
    // (The Mathworks, 2012)). This data is generated using the DOPRI56 numerical integrator.
    std::string pathToDiscreteEventBenchmarkDatafile = tudat::input_output::getTudatRootPath( )
            + "/Mathematics/NumericalIntegrators/UnitTests"
            + "/matlabOutputDiscreteEventRungeKutta56DormandPrince.txt";

    // Run Matlab integrator tests.
    matlab_numerical_integrator_tests::runMatlabNumericalIntegratorTests(
                pathToBenchmarkDatafile, 1.0e-13, 1.0e-10,
                RungeKuttaCoefficients::get( RungeKuttaCoefficients::rungeKuttaFehlberg56 ),
                pathToDiscreteEventBenchmarkDatafile );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

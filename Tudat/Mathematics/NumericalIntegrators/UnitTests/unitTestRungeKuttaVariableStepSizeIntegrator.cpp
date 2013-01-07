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
 *      120203    B. Tong Minh      Copied RungeKutta4Stepsize unit test.
 *      120207    K. Kumar          Adapted to use modified benchmark functions in Tudat Core.
 *      120213    K. Kumar          Modified getCurrentInterval() to getIndependentVariable();
 *                                  transferred to Boost unit test framework.
 *      120321    K. Kumar          Updated (Burden and Faires, 2011) benchmark function call.
 *      120328    K. Kumar          Removed specific integrator tests; added compiler tests; added
 *                                  test of runtime error when minimm step size is exceeded.
 *      120331    B. Tong Minh      Added typeid check for RungeKuttaVariableStepSizeIntegratorXd
 *                                  typedef; modified minimum step size exceeded unit tests to use
 *                                  custom exception object.
 *
 *    References
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole, 2001.
 *
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <typeinfo>

#include <boost/exception/all.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h"
#include "Tudat/Mathematics/NumericalIntegrators/UnitTests/numericalIntegratorTestFunctionSuite.h"

namespace tudat
{
namespace unit_tests
{

using tudat::numerical_integrators::RungeKuttaCoefficients;
using tudat::numerical_integrators::RungeKuttaVariableStepSizeIntegrator;
using tudat::numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd;
using numerical_integrator_test_function_suite::computeZeroStateDerivative;

BOOST_AUTO_TEST_SUITE( test_runge_kutta_variable_step_size_integrator )

//! Test different types of states and state derivatives.
BOOST_AUTO_TEST_CASE( testCompilerErrors )
{
    // Case 1: test the VectorXd typdef. There is no need to explicitly test anything, since we're
    //         looking for compiler errors.
    {
        RungeKuttaVariableStepSizeIntegratorXd integrator(
                    RungeKuttaCoefficients( ),
                    &computeZeroStateDerivative,
                    0.0,
                    Eigen::Vector3d::Zero( ),
                    0.01,
                    std::numeric_limits< double >::infinity( ),
                    10.0 * std::numeric_limits< double >::epsilon( ),
                    10.0 * std::numeric_limits< double >::epsilon( ) );

        // Test integrateTo() function.
        integrator.integrateTo( 10.0, 0.1 );

        // Test performIntegrationStep() function.
        integrator.performIntegrationStep( 0.1 );

        // Test rollbackToPreviousState() function.
        integrator.rollbackToPreviousState( );

        // Check if the default types are correct.
        BOOST_CHECK_EQUAL(
                    std::string( typeid( integrator.getCurrentIndependentVariable( ) ).name( ) ),
                    std::string( typeid( double ).name( ) ) );
        BOOST_CHECK_EQUAL( std::string( typeid( integrator.getCurrentState( ) ).name( ) ),
                           std::string( typeid( Eigen::VectorXd ).name( ) ) );
    }

    // Case 2: test different types of states and state derivatives. There is no need to explicitly
    //         test anything, since we're looking for compiler errors.
    {
        RungeKuttaVariableStepSizeIntegrator < double, Eigen::Vector3d, Eigen::VectorXd >
                integrator( RungeKuttaCoefficients( ),
                            &computeZeroStateDerivative,
                            0.0,
                            Eigen::Vector3d::Zero( ),
                            0.01,
                            std::numeric_limits< double >::infinity( ),
                            10.0 * std::numeric_limits< double >::epsilon( ),
                            10.0 * std::numeric_limits< double >::epsilon( ) );

        // Test integrateTo() function.
        integrator.integrateTo( 10.0, 0.1 );

        // Test performIntegrationStep() function.
        integrator.performIntegrationStep( 0.1 );

        // Test rollbackToPreviousState() function.
        integrator.rollbackToPreviousState( );
    }

    // Case 3: test same types of states and state derivatives (MatrixXd). There is no need to
    //         explicitly test anything, since we're looking for compiler errors.
    {
        RungeKuttaVariableStepSizeIntegrator < double, Eigen::MatrixXd, Eigen::MatrixXd >
                integrator( RungeKuttaCoefficients( ),
                            &computeZeroStateDerivative,
                            0.0,
                            Eigen::Vector3d::Zero( ),
                            0.01,
                            std::numeric_limits< double >::infinity( ),
                            10.0 * std::numeric_limits< double >::epsilon( ),
                            10.0 * std::numeric_limits< double >::epsilon( ) );

        // Test integrateTo() function.
        integrator.integrateTo( 10.0, 0.1 );

        // Test performIntegrationStep() function.
        integrator.performIntegrationStep( 0.1 );

        // Test rollbackToPreviousState() function.
        integrator.rollbackToPreviousState( );
    }
}

//! Test that exceeding minimum step size throws a runtime error.
BOOST_AUTO_TEST_CASE( testMinimumStepSizeRuntimeError )
{
    RungeKuttaVariableStepSizeIntegratorXd integrator(
                RungeKuttaCoefficients( ),
                &computeZeroStateDerivative,
                0.0,
                Eigen::Vector3d::Zero( ),
                100.0,
                std::numeric_limits< double >::infinity( ),
                std::numeric_limits< double >::epsilon( ),
                std::numeric_limits< double >::epsilon( ) );

    // Case 1: test that minimum step size is exceeded when using integrateTo().
    {
        // Declare boolean flag to test if minimum step size is exceed for integrateTo().
        bool isMinimumStepSizeExceededForIntegrateTo = false;

        // Try integrateTo(), which should result in a runtime error.
        try
        {
            // Test integrateTo() function.
            integrator.integrateTo( 10.0, 0.1 );
        }

        // Catch the expected runtime error, and set the boolean flag to true.
        catch ( RungeKuttaVariableStepSizeIntegratorXd::MinimumStepSizeExceededError
                minimumStepSizeExceededError )
        {
            isMinimumStepSizeExceededForIntegrateTo = true;
            BOOST_CHECK_EQUAL( minimumStepSizeExceededError.minimumStepSize, 100.0 );
        }

        // Check that the minimum step size was indeed exceeded.
        BOOST_CHECK( isMinimumStepSizeExceededForIntegrateTo );
    }

    // Case 2: test that minimum step size is exceeded when using performIntegrationStep().
    {
        // Declare boolean flag to test if minimum step size is exceed for
        // performIntegrationStep().
        bool isMinimumStepSizeExceededForPerformIntegrationStep = false;

        // Try performIntegrationStep(), which should result in a runtime error.
        try
        {
            // Test performIntegrationStep() function.
            integrator.performIntegrationStep( 0.1 );
        }

        // Catch the expected runtime error, and set the boolean flag to true.
        catch ( RungeKuttaVariableStepSizeIntegratorXd::MinimumStepSizeExceededError
                minimumStepSizeExceededError )
        {
            isMinimumStepSizeExceededForPerformIntegrationStep = true;
            BOOST_CHECK_EQUAL( minimumStepSizeExceededError.minimumStepSize, 100.0 );
        }

        // Check that the minimum step size was indeed exceeded.
        BOOST_CHECK( isMinimumStepSizeExceededForPerformIntegrationStep );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole, 2001.
 *
 *    Notes
 *      This file doesn't test any specific Runge-Kutta-type integrators, but rather some general
 *      functionality adopted in the the RungeKuttaVariableStepSizeIntegrator class, applicable to
 *      all Runge-Kutta-type integrators.
 *
 *      It should be noted that the getCurrentStateDerivatives() member function is not tested in a
 *      fully generic manner at the moment; the test is setup specifically based on the 
 *      Runge-Kutta-Fehlberg 4(5) (RKF45) integrator. A more comprehensive test should be designed
 *      to ensure that the member function performs as desired regardless of the coefficient set
 *      chosen.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN


#include <boost/test/unit_test.hpp>

#include "tudat/math/integrators/rungeKuttaVariableStepSizeIntegrator.h"
#include "tudat/math/integrators/rungeKutta4Integrator.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/basics/testMacros.h"
#include "tudat/math/integrators/numericalIntegratorTestFunctions.h"

#include <limits>
#include <string>
#include <typeinfo>

namespace tudat
{
namespace unit_tests
{

using numerical_integrators::RungeKuttaCoefficients;
using numerical_integrators::RungeKuttaVariableStepSizeIntegrator;
using numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd;
using numerical_integrator_test_functions::computeZeroStateDerivative;

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

//! Test if the state derivative evaliations are properly returned.
BOOST_AUTO_TEST_CASE( testStateDerivativeRetrievalFunction )
{
    using namespace numerical_integrators;
    using namespace unit_tests::numerical_integrator_test_functions;

    // This test is based on the Runge-Kutta-Fehlberg 4(5) coefficient set, hence the test does not
    // robustly ensure that the getCurrentStateDerivatives() works correctly for any given 
    // coefficient set currently. This test is more of an preliminary check that the function 
    // performs as required, with the extrapolation that it is likely to perform consistently in
    // this manner for any Runge-Kutta-type coefficient set.

    // Create test integrator.
    RungeKuttaVariableStepSizeIntegratorXd integrator(
                RungeKuttaCoefficients::get( CoefficientSets::rungeKuttaFehlberg45 ),
                &computeVanDerPolStateDerivative,
                0.0,
                ( Eigen::VectorXd( 2 ) << 1.0, 2.0 ).finished( ),
                0.0, 10.0, 1.0E-8, 1.0E-8 );

    // Perform first integration step.
    integrator.performIntegrationStep( 1.0 );

    // Retrieve time and state after first integration step.
    const double previousTime = integrator.getCurrentIndependentVariable( );
    const Eigen::VectorXd previousState = integrator.getCurrentState( );

    // Perform additional integration step to verify that there are no problems with
    // getCurrentStateDerivatives after >1 iterations (i.e not being reset etc.).
    integrator.performIntegrationStep( integrator.getNextStepSize( ) );

    // Get current time and previous time step.
    const double currentTime = integrator.getCurrentIndependentVariable( );
    const double stepSize = currentTime - previousTime;

    // Retrieve state derivative values used in previous time step.
    std::vector< Eigen::VectorXd > stateDerivatives = integrator.getCurrentStateDerivatives( );

    // Check size of state derivative vector. This test is specifically set up to test the for the 
    // number of stages in the RKF45 integrator.
    BOOST_CHECK_EQUAL( stateDerivatives.size( ), 6 );

    // Perform manual state derivative evaluations at previous time step using RKF45 method
    RungeKuttaCoefficients rkf45Coefficients = RungeKuttaCoefficients::get(
                rungeKuttaFehlberg45 );
    std::vector< Eigen::VectorXd > directStateDerivativeValues;

    for ( int stage = 0; stage < rkf45Coefficients.cCoefficients.rows( ); stage++ )
    {
        // Compute the intermediate state to pass to the state derivative for this stage.
        Eigen::VectorXd intermediateState( previousState );

        // Compute the intermediate state.
        for ( int column = 0; column < stage; column++ )
        {
            intermediateState += stepSize * rkf45Coefficients.aCoefficients( stage, column )
                    * directStateDerivativeValues[ column ];
        }

        // Compute state derivative.
        directStateDerivativeValues.push_back(
                    computeVanDerPolStateDerivative(
                        previousTime +
                        rkf45Coefficients.cCoefficients( stage ) * stepSize,
                        intermediateState ) );

        // Check if manual result matched result from NumericalIntegrator.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( directStateDerivativeValues.at( stage ),
                                           stateDerivatives.at( stage ),
                                           1.0E-15 );
    }
}

//! Test if integtrateTo function works for variable step size where last step is modified.
BOOST_AUTO_TEST_CASE( testVariableStepIntegrateToFunction )
{
    using namespace numerical_integrators;
    using namespace unit_tests::numerical_integrator_test_functions;

    // In this test, the integrateTo function is used with variable step size integrator, where the
    // step size is modified by the last step size, to ensure the correct operation of the
    // integrator in this case.

    // Create variable step size integrator.
    RungeKuttaVariableStepSizeIntegratorXd integrator(
                RungeKuttaCoefficients::get( CoefficientSets::rungeKuttaFehlberg45 ),
                &computeSecondNonAutonomousModelStateDerivative,
                0.0,
                ( Eigen::VectorXd( 1 ) << 1.0 ).finished( ),
                0.0, 10.0, 1.0E-10, 1.0E-10 );

    // Use integrateTo function
    Eigen::VectorXd integratedValue = integrator.integrateTo( 0.5, 1.0 );
    double currentTime = integrator.getCurrentIndependentVariable( );

    // Check if current time is correct.
    BOOST_CHECK_CLOSE_FRACTION( currentTime, 0.5, std::numeric_limits< double >::epsilon( ) );

    // Create integrator to see if performing a single step with the given settings will indeed
    // result in the step size being adapted.
    RungeKuttaVariableStepSizeIntegratorXd verificationIntegrator(
                RungeKuttaCoefficients::get( CoefficientSets::rungeKuttaFehlberg45 ),
                &computeSecondNonAutonomousModelStateDerivative,
                0.0,
                ( Eigen::VectorXd( 1 ) << 1.0 ).finished( ),
                0.0, 10.0, 1.0E-10, 1.0E-10 );

    // Check if single step size of 0.5 will be adapted.
    verificationIntegrator.performIntegrationStep( 0.5 );
    BOOST_CHECK_EQUAL( ( verificationIntegrator.getCurrentIndependentVariable( ) < 0.5 *
                         ( 1.0 - 10.0 * std::numeric_limits< double >::epsilon( ) ) ), true );

    // Use a fixed step size integrator to check the original result of integrateTo
    RungeKutta4IntegratorXd fixedStepSizeIntegrator(
                &computeSecondNonAutonomousModelStateDerivative,0.0,
                ( Eigen::VectorXd( 1 ) << 1.0 ).finished( ) );
    Eigen::VectorXd fixedStepIntegratedValue = integrator.integrateTo( 0.5, 0.01 );
    BOOST_CHECK_CLOSE_FRACTION( fixedStepIntegratedValue.x( ), integratedValue.x( ), 1.0E-10 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

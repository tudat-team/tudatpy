/*    Copyright (c) 2010-2017, Delft University of Technology
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
 *      The MathWorks, Inc. Symbolic Math Toolbox, 2012.
 *
 */

#define BOOST_TEST_MAIN

#include <iostream>

#include <limits>

#include <boost/test/unit_test.hpp>

#include <Tudat/Mathematics/NumericalIntegrators/UnitTests/numericalIntegratorTestFunctions.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Mathematics/NumericalIntegrators/bulirschStoerVariableStepsizeIntegrator.h"

namespace tudat
{
namespace unit_tests
{

const double TUDAT_MACHINE_PRECISION = std::numeric_limits< double >::epsilon( );

BOOST_AUTO_TEST_SUITE( test_bulirsch_stoer_variable_step_size_integrator )

using namespace tudat::numerical_integrators;

//! Test the validity of the BurlischStoerVariableStepSize integrator.
/*!
 * Tests the validity of the RungeKuttaVariableStepSize integrator.
 * \param sequence The rational function sequence used for the Bulirsch-Stoer integrator being
 *          tested.
 * \param stateDerivativeFunction Function pointer to the state derivative function.
 * \param intervalStart The start of the integration interval.
 * \param intervalEnd The end of the integration interval.
 * \param stepSize The step size to take.
 * \param initialState The initial state.
 * \param expectedState Expected final state.
 * \param tolerance Tolerance when comparing.
 * \return True if actual final state equals the expected final state within the specified
 *          tolerance.
 */
bool testValidityOfBulirschStoerVariableStepSizeIntegrator(
        const std::vector< unsigned int > sequence,
        const BulirschStoerVariableStepSizeIntegratorXd::StateDerivativeFunction&
        stateDerivativeFunction,
        const double intervalStart, const double intervalEnd, const double stepSize,
        const Eigen::VectorXd& initialState, const Eigen::VectorXd expectedState,
        const double tolerance )
{
    // Create forward BulirschStoerVariableStepSize integrator.
    {
        std::cout<<"T1"<<std::endl;
        BulirschStoerVariableStepSizeIntegratorXd integrator( sequence, stateDerivativeFunction,
                                                              intervalStart, initialState,
                                                              1.0e-15, tolerance / 10.0 );
        std::cout<<"T1"<<std::endl;


        Eigen::VectorXd finalState = integrator.integrateTo( intervalEnd, stepSize );

        // Compute differences between computed and expected interval end and generate
        // cerr statement if test fails.
        if ( std::fabs( integrator.getCurrentIndependentVariable( ) - intervalEnd ) / intervalEnd >
             10.0 * TUDAT_MACHINE_PRECISION )
        {
            return false;
        }

        // Compute differences between computed and expected results and generate
        // cerr statement if test fails.
        if ( !expectedState.isApprox( finalState, tolerance ) )
        {
            return false;
        }
    }
    std::cout<<"T1"<<std::endl;


    // Try the same again, but in two steps.
    {
        std::cout<<"T2"<<std::endl;

        BulirschStoerVariableStepSizeIntegratorXd integrator( sequence, stateDerivativeFunction,
                                                              intervalStart, initialState,
                                                              1.0e-15, tolerance / 10.0 );

        const double intermediateIndependentVariable
                = intervalStart + ( intervalEnd - intervalStart ) / 2.0;

        const Eigen::VectorXd intermediateState = integrator.integrateTo(
                    intermediateIndependentVariable, stepSize );

        // Compute differences between computed and expected interval end.
        if ( std::fabs( integrator.getCurrentIndependentVariable( )
                       - intermediateIndependentVariable ) /
             intermediateIndependentVariable > TUDAT_MACHINE_PRECISION )
        {
            return false;
        }

        std::cout<<"T2"<<std::endl;

        // Integrate to the end
        Eigen::VectorXd finalState = integrator.integrateTo( intervalEnd, stepSize );

        // Compute differences between computed and expected interval end.
        if ( std::fabs( integrator.getCurrentIndependentVariable( ) - intervalEnd ) / intervalEnd >
             10.0 * TUDAT_MACHINE_PRECISION )
        {
            return false;
        }

        // Compute differences between computed and expected results and generate
        // cerr statement if test fails.
        if ( !expectedState.isApprox( finalState, tolerance ) )
        {
            return false;
        }

        integrator.performIntegrationStep( stepSize );
        if ( !integrator.rollbackToPreviousState( ) )
        {
            return false;
        }

        // No need to check machine precision, because this interval is stored exact.
        if ( std::fabs( integrator.getCurrentIndependentVariable( ) - intervalEnd ) / intervalEnd >
             10.0 * TUDAT_MACHINE_PRECISION )
        {
            return false;
        }

        // This result should be exactly the same.
        if ( integrator.getCurrentState( ) != finalState )
        {
            return false;
        }

        if ( integrator.rollbackToPreviousState( ) )
        {
            return false;
        }
    }

    return true;
}

//! Test different types of states and state derivatives.
/*!
 * Tests if different types of states and state derivatives work. If something
 * is broken, then a compile time error will be generated.
 * \return Unconditionally true
 */
bool testDifferentStateAndStateDerivativeTypes( )
{
    using tudat::unit_tests::numerical_integrator_test_functions::computeZeroStateDerivative;
    tudat::numerical_integrators::BulirschStoerVariableStepSizeIntegrator
            < double, Eigen::Vector3d, Eigen::VectorXd > integrator(
                getBulirschStoerStepSequence( ),  &computeZeroStateDerivative,
                0.0, Eigen::Vector3d::Zero( ) );

    integrator.integrateTo( 0.0, 0.1 );

    // No need to test anything, this is just to check compile time errors.
    return true;
}

//! Test the validity of the BulirschStoerVariableStepSize integrator.
/*!
 * Tests the validity of the BulirschStoerVariableStepSize integrator.
 * \param sequence Bulirsch-Stoer, variable-step size rational function sequence.
 * \return True if actual final state equals the expected final state within the specified
 *          tolerance.
 */
bool testValidityOfBulirschStoerVariableStepSizeIntegrator(
        const std::vector< unsigned int >& sequence )
{
    // Test result initialised to false.
    bool testBulirschStoerVariableStepSizeIsOk = true;

    using namespace tudat::unit_tests;
//    std::map< BenchmarkFunctions, BenchmarkFunction >& benchmarkFunctions =
//            getBenchmarkFunctions( );

//    // Case 1: test with x_dot = 0, which results in x_f = x_0.
//    {
//        testBulirschStoerVariableStepSizeIsOk
//                &= testValidityOfBulirschStoerVariableStepSizeIntegrator(
//                    sequence,
//                    &unit_tests::numerical_integrator_test_functions::computeZeroStateDerivative,
//                    0.0,
//                    100.0, 0.2,
//                    Eigen::Vector3d::Zero( ),
//                    Eigen::Vector3d::Zero( ),
//                    TUDAT_MACHINE_PRECISION );
//    }

//    // Case 2: test with x_dot = 1, which results in x_f = x_0 + t_f.
//    {
//        std::cout<<"test 2"<<std::endl;
//        testBulirschStoerVariableStepSizeIsOk
//                &= testValidityOfBulirschStoerVariableStepSizeIntegrator(
//                    sequence,
//                    &unit_tests::numerical_integrator_test_functions::computeConstantStateDerivative,
//                    0.0,
//                    100.0, 0.2,
//                    Eigen::Vector3d::UnitX( ),
//                    Eigen::Vector3d::UnitX( ) * 101.0,
//                    1000.0 * TUDAT_MACHINE_PRECISION );
//        std::cout<<"test 2 done"<<std::endl;

//    }

    // Case 3: test with x_dot = x, which results in x_f = x0 * exp( t_f )
    {
        std::cout<<"test 3"<<std::endl;

        testBulirschStoerVariableStepSizeIsOk
                &= testValidityOfBulirschStoerVariableStepSizeIntegrator(
                    sequence,
                    &unit_tests::numerical_integrator_test_functions::computeVanDerPolStateDerivative,
                    0.0,
                    10.0, 1.0,
                    Eigen::Vector2d::UnitX( ),
                    Eigen::Vector2d::UnitX( ) * std::exp( 10.0 ), 1.0e-14 );
        std::cout<<"test 3 done"<<std::endl;

    }

//    // Case 4: test with an example from Burden and Faires.
//    {
//        testBulirschStoerVariableStepSizeIsOk
//                &= testValidityOfBulirschStoerVariableStepSizeIntegrator(
//                    sequence,
//                    benchmarkFunctions[ BurdenAndFaires ]
//                    .pointerToStateDerivativeFunction_,
//                    benchmarkFunctions[ BurdenAndFaires ].intervalStart_,
//                    benchmarkFunctions[ BurdenAndFaires ].intervalEnd_, 0.25,
//                    benchmarkFunctions[ BurdenAndFaires ].initialState_,
//                    benchmarkFunctions[ BurdenAndFaires ].finalState_, 1.0e-11 );
//    }

    return testBulirschStoerVariableStepSizeIsOk;
}

BOOST_AUTO_TEST_CASE( testBulirschStoerVariableStepSizeIntegrator )
{
    // Case 1: test if difference in type between state and state derivative works.
    //BOOST_CHECK( testDifferentStateAndStateDerivativeTypes( ) );

    // Case 2: test if Bulirsch-Stoer sequence works.
    BOOST_CHECK( testValidityOfBulirschStoerVariableStepSizeIntegrator(
                     getBulirschStoerStepSequence( ) ) );

//    // Case 2: test if Deufelhard sequence works.
//    BOOST_CHECK( testValidityOfBulirschStoerVariableStepSizeIntegrator(
//                     BulirschStoerRationalFunctionSequences::get(
//                         BulirschStoerRationalFunctionSequences::deufelhard ) ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

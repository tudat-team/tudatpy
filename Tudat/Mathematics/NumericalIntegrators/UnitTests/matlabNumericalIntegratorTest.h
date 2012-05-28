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
 *      120329    K. Kumar          File created.
 *      120404    K. Kumar          Added unit test for discrete events.
 *
 *    References
 *      The Mathworks, Inc. RKF78, Symbolic Math Toolbox, 2012.
 *
 */

#ifndef TUDAT_MATLAB_NUMERICAL_INTEGRATOR_TEST_H
#define TUDAT_MATLAB_NUMERICAL_INTEGRATOR_TEST_H

#include <string>

#include <boost/bind.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/matrixTextFileReader.h"
#include "Tudat/Mathematics/NumericalIntegrators/UnitTests/numericalIntegratorTestFunctionSuite.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h"

namespace tudat
{
namespace unit_tests
{
namespace matlab_numerical_integrator_tests
{

using tudat::mathematics::numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd;
using tudat::mathematics::numerical_integrators::RungeKuttaCoefficients;

void runMatlabNumericalIntegratorTests( const std::string& pathToBenchmarkDataFile,
                                        const double singleStepTestTolerance,
                                        const double fullIntegrationTestTolerance,
                                        const RungeKuttaCoefficients& coefficients,
                                        const std::string& pathToDiscreteEventBenchmarkDataFile = "" )
{
    // Store benchmark data in matrix.
    const Eigen::MatrixXd matlabBenchmarkData =
            input_output::readMatrixFromFile( pathToBenchmarkDataFile, " ," );

    // Store discrete event benchmark data in matrix.
    const Eigen::MatrixXd matlabDiscreteEventBenchmarkData =
            input_output::readMatrixFromFile( pathToDiscreteEventBenchmarkDataFile, " ," );

    // Set parameters of integration.
    // This should to be added to the benchmark data file and parsed accordingly once the Tudat
    // parser architecture has been added to the library. The relative and absolute error tolerance
    // are set large such that the condition in the numerical integrator is always satisfied
    // to accept the given step size for the next step.
    const double initialTime = -1.0;
    const double finalTime = 1.0;
    const Eigen::VectorXd initialState = Eigen::VectorXd::Constant( 1, 0.5 );
    const double minimumStepSize = 1.0e-16;
    const double maximumStepSize = 10.0;
    const double relativeErrorTolerance = 1.0e-15;
    const double absoluteErrorTolerance = 10.0 * 1.0e-16;

    // Set index of last row of matrix containing benchmark data.
    const int indexLastRowBenchmarkData = matlabBenchmarkData.rows( ) - 1;

    // Set up tolerances so that the step size provided is always accepted by the integrator for
    // the test cases that need that behavior.
    const double infiniteRelativeErrorTolerance = std::numeric_limits< double >::infinity( );
    const double infiniteAbsoluteErrorTolerance = std::numeric_limits< double >::infinity( );

    // Set discrete event times.
    const double TIME_OF_DISCRETE_EVENT_1 = -0.5;
    const double TIME_OF_DISCRETE_EVENT_2 = 0.0;
    const double TIME_OF_DISCRETE_EVENT_3 = 0.5;

    // Set column indices in matrices.
    const int TIME_COLUMN_INDEX = 0;
    const int STATE_COLUMN_INDEX = 1;
    const int FIRST_ROW = 0;
    const int SECOND_ROW = 1;

    using numerical_integrator_test_function_suite::computeNonAutonomousModelStateDerivative;

    // Case 1: Use integrateTo() to integrate one step forward in time and check results against
    //         benchmark data from (The Mathworks, 2012).
    {
        // Declare integrator with all necessary settings.
        RungeKuttaVariableStepSizeIntegratorXd integrator(
                    coefficients,
                    &computeNonAutonomousModelStateDerivative,
                    initialTime,
                    initialState,
                    minimumStepSize,
                    maximumStepSize,
                    infiniteRelativeErrorTolerance,
                    infiniteAbsoluteErrorTolerance );

        // Check that it is now not possible to roll back.
        BOOST_CHECK( !integrator.rollbackToPreviousState( ) );

        // Use integrateTo() to integrate one step.
        integrator.integrateTo( matlabBenchmarkData( SECOND_ROW, TIME_COLUMN_INDEX ),
                                matlabBenchmarkData( SECOND_ROW, TIME_COLUMN_INDEX )
                                - matlabBenchmarkData( FIRST_ROW, TIME_COLUMN_INDEX ) );

        // Check if the expected time at the end of the integration step matches the required
        // time. This test should be exact.
        BOOST_CHECK_EQUAL( matlabBenchmarkData( SECOND_ROW, TIME_COLUMN_INDEX ),
                           integrator.getCurrentIndependentVariable( ) );

        // Check if expected state at end of integration step matches required state.
        BOOST_CHECK_CLOSE_FRACTION( matlabBenchmarkData( SECOND_ROW, STATE_COLUMN_INDEX ),
                                    integrator.getCurrentState( )( 0 ), singleStepTestTolerance );

        // Roll back to the previous step. This should be possible since the integrateTo() function
        // was called above.
        BOOST_CHECK( integrator.rollbackToPreviousState( ) );

        // Check that the rolled back time is as required. This test should be exact.
        BOOST_CHECK_EQUAL( matlabBenchmarkData( FIRST_ROW, TIME_COLUMN_INDEX ),
                           integrator.getCurrentIndependentVariable( ) );

        // Check that the rolled back state is as required. This test should be exact.
        BOOST_CHECK_EQUAL( matlabBenchmarkData( FIRST_ROW, STATE_COLUMN_INDEX ),
                           integrator.getCurrentState( )( 0 ) );

        // Check that it is now not possible to roll back.
        BOOST_CHECK( !integrator.rollbackToPreviousState( ) );
    }

    // Case 2: Use performIntegrationstep() to integrate to final time in multiple steps and check
    //         results against benchmark data from (The Mathworks, 2012). Unfortunately, the
    //         adaptive step size method using in Matlab is not well-documented, so it circumvent
    //         this, the Tudat numerical integrator tested here is provided with the steps that
    //         must be taken, based on the time-series data in the data file. The integrator is
    //         set up such that the step given based on this time-series data is always accepted.
    {
        RungeKuttaVariableStepSizeIntegratorXd integrator(
                    coefficients,
                    &computeNonAutonomousModelStateDerivative,
                    initialTime,
                    initialState,
                    minimumStepSize,
                    maximumStepSize,
                    infiniteRelativeErrorTolerance,
                    infiniteAbsoluteErrorTolerance );

        // Declare last time and state.
        double lastTime = TUDAT_NAN;
        Eigen::VectorXd lastState = Eigen::VectorXd::Zero( 1 );

        // Loop through all the integration steps in the benchmark data.
        for ( int i = 1; i < matlabBenchmarkData.rows( ); i++ )
        {
            // Set the step size based on the output data generated with Matlab.
            double stepSize = matlabBenchmarkData( i, TIME_COLUMN_INDEX )
                    - matlabBenchmarkData( i - 1, TIME_COLUMN_INDEX );

            // Store last time and state.
            lastTime = integrator.getCurrentIndependentVariable( );
            lastState = integrator.getCurrentState( );

            // Perform the integration step.
            integrator.performIntegrationStep( stepSize );

            // Check if the expected time at the end of the integration step matches the required
            // time. This test should be accurate to effectively machine precision.
            BOOST_CHECK_CLOSE_FRACTION( matlabBenchmarkData( i, TIME_COLUMN_INDEX ),
                                        integrator.getCurrentIndependentVariable( ),
                                        10.0 * std::numeric_limits< double >::epsilon( ) );

            // Check if expected state at end of integration step matches required state.
            // The reason this test has a different, higher tolerance than for a single integration
            // step is because the acummulated error builds up, since this test is based off of
            // Matlab time-series. The individual steps are accurate to the tolerance used in
            // Case 1.
            BOOST_CHECK_CLOSE_FRACTION( matlabBenchmarkData( i, STATE_COLUMN_INDEX ),
                                        integrator.getCurrentState( )( 0 ),
                                        fullIntegrationTestTolerance );
        }

        // Roll back to the previous step. This should be possible since the
        // performIntegrationStep() function was called above.
        BOOST_CHECK( integrator.rollbackToPreviousState( ) );

        // Check that the rolled back time is as required. This test should be exact.
        BOOST_CHECK_EQUAL( lastTime, integrator.getCurrentIndependentVariable( ) );

        // Check that the rolled back state is as required. This test should be exact.
        BOOST_CHECK_EQUAL( lastState( 0 ), integrator.getCurrentState( )( 0 ) );

        // Check that it is now not possible to roll back.
        BOOST_CHECK( !integrator.rollbackToPreviousState( ) );
    }

    // Case 3: Use performIntegrationstep() to integrate to initial time in multiple (backwards)
    //         steps and check results against benchmark data from (The Mathworks, 2012).
    //         Unfortunately, the adaptive step size method using in Matlab is not well-documented,
    //         so it circumvent this, the Tudat numerical integrator tested here is provided with
    //         the steps that must be taken, based on the time-series data in the data file. The
    //         integrator is set up such that the step given based on this time-series data is
    //         always accepted.
    {
        RungeKuttaVariableStepSizeIntegratorXd integrator(
                    coefficients,
                    &computeNonAutonomousModelStateDerivative,
                    matlabBenchmarkData( indexLastRowBenchmarkData, TIME_COLUMN_INDEX ),
                    matlabBenchmarkData.block( indexLastRowBenchmarkData,
                                               STATE_COLUMN_INDEX, 1, 1 ),
                    minimumStepSize,
                    maximumStepSize,
                    infiniteRelativeErrorTolerance,
                    infiniteAbsoluteErrorTolerance );

        // Declare last time and state.
        double lastTime = TUDAT_NAN;
        Eigen::VectorXd lastState = Eigen::VectorXd::Zero( 1 );

        for ( int i = indexLastRowBenchmarkData; i > 0; --i )
        {
            // Set the step size based on the output data generated with Matlab.
            double stepSize = matlabBenchmarkData( i - 1, TIME_COLUMN_INDEX )
                    - matlabBenchmarkData( i, TIME_COLUMN_INDEX );

            // Store last time and state.
            lastTime = integrator.getCurrentIndependentVariable( );
            lastState = integrator.getCurrentState( );

            // Perform the integration step.
            integrator.performIntegrationStep( stepSize );

            // Check if the expected time at the end of the integration step matches the required
            // time. This test should be accurate to effectively machine precision.
            BOOST_CHECK_CLOSE_FRACTION( matlabBenchmarkData( i - 1, TIME_COLUMN_INDEX ),
                                        integrator.getCurrentIndependentVariable( ),
                                        10.0 * std::numeric_limits< double >::epsilon( ) );

            // Check if expected state at end of integration step matches required state.
            // The reason this test has a different, higher tolerance than for a single integration
            // step is because the acummulated error builds up, since this test is based off of
            // Matlab time-series. The individual steps are accurate to the tolerance used in
            // Case 1.
            BOOST_CHECK_CLOSE_FRACTION( matlabBenchmarkData( i - 1, STATE_COLUMN_INDEX ),
                                        integrator.getCurrentState( )( 0 ),
                                        fullIntegrationTestTolerance );
        }

        // Check that it is possible to roll back to the previous step.
        BOOST_CHECK( integrator.rollbackToPreviousState( ) );

        // Check that the rolled back time is as required. This test should be exact.
        BOOST_CHECK_EQUAL( lastTime, integrator.getCurrentIndependentVariable( ) );

        // Check that the rolled back state is as required. This test should be exact.
        BOOST_CHECK_EQUAL( lastState( 0 ), integrator.getCurrentState( )( 0 ) );

        // Check that it is now not possible to roll back to the previous step.
        BOOST_CHECK( !integrator.rollbackToPreviousState( ) );
    }

    // Case 4: Use integrateTo() to integrate to final time and and check results against benchmark
    //         data from (The Mathworks, 2012).
    {
        // Declare integrator with all necessary settings.
        RungeKuttaVariableStepSizeIntegratorXd integrator(
                    coefficients,
                    &computeNonAutonomousModelStateDerivative,
                    initialTime,
                    initialState,
                    minimumStepSize,
                    maximumStepSize,
                    relativeErrorTolerance,
                    absoluteErrorTolerance );

        // Use integrateTo() to integrate to final time.
        integrator.integrateTo( finalTime,
                                matlabBenchmarkData( SECOND_ROW, TIME_COLUMN_INDEX )
                                - matlabBenchmarkData( FIRST_ROW, TIME_COLUMN_INDEX ) );

        // Check if the expected time at the end of the integration step matches the required
        // time. This test should be exact.
        BOOST_CHECK_EQUAL( matlabBenchmarkData( indexLastRowBenchmarkData, TIME_COLUMN_INDEX ),
                           integrator.getCurrentIndependentVariable( ) );

        // Check if expected state at end of integration step matches required state.
        // The reason this test has a different, higher tolerance than for a single integration
        // step is because the acummulated error builds up, since this test is based off of
        // Matlab time-series. The individual steps are accurate to the tolerance used in
        // Case 1.
        BOOST_CHECK_CLOSE_FRACTION( matlabBenchmarkData( indexLastRowBenchmarkData,
                                                         STATE_COLUMN_INDEX ),
                                    integrator.getCurrentState( )( 0 ),
                                    fullIntegrationTestTolerance );

        // Roll back to the previous step. This should be possible since the integrateTo() function
        // was called above.
        BOOST_CHECK( integrator.rollbackToPreviousState( ) );

        // Check that it is now not possible to roll back.
        BOOST_CHECK( !integrator.rollbackToPreviousState( ) );
    }

    // Case 5: Use performIntegrationstep() to integrate to final time in multiple steps, including
    //         discrete events and check results against benchmark data from (The Mathworks, 2012).
    //         Unfortunately, the adaptive step size method using in Matlab is not well-documented,
    //         so it circumvent this, the Tudat numerical integrator tested here is provided with
    //         the steps that must be taken, based on the time-series data in the data files. The
    //         integrator is set up such that the step given based on this time-series data is
    //         always accepted. The discrete events are enabled by modifying the state
    //         instantaneously. These discrete events cannot be rolled back.
    {
        RungeKuttaVariableStepSizeIntegratorXd integrator(
                    coefficients,
                    &computeNonAutonomousModelStateDerivative,
                    initialTime,
                    initialState,
                    minimumStepSize,
                    maximumStepSize,
                    infiniteRelativeErrorTolerance,
                    infiniteAbsoluteErrorTolerance );

        // Declare last time and state.
        double lastTime = TUDAT_NAN;
        Eigen::VectorXd lastState = Eigen::VectorXd::Zero( 1 );

        for ( int i = 1; i < matlabDiscreteEventBenchmarkData.rows( ); i++ )
        {
            // Set the step size based on the output data generated with Matlab.
            double stepSize = matlabDiscreteEventBenchmarkData( i, TIME_COLUMN_INDEX )
                    - matlabDiscreteEventBenchmarkData( i - 1, TIME_COLUMN_INDEX );

            // Store last time and state.
            lastTime = integrator.getCurrentIndependentVariable( );
            lastState = integrator.getCurrentState( );

            // Perform the integration step.
            integrator.performIntegrationStep( stepSize );

            // Check if the expected time at the end of the integration step matches the required
            // time. This test should be accurate to effectively machine precision.
            if ( std::fabs( matlabDiscreteEventBenchmarkData( i, TIME_COLUMN_INDEX ) )
                 < std::numeric_limits< double >::epsilon( ) )
            {
                BOOST_CHECK_SMALL( integrator.getCurrentIndependentVariable( ),
                                   std::numeric_limits< double >::epsilon( ) );
            }
            else
            {
                BOOST_CHECK_CLOSE_FRACTION(
                            matlabDiscreteEventBenchmarkData( i, TIME_COLUMN_INDEX ),
                            integrator.getCurrentIndependentVariable( ),
                            10.0 * std::numeric_limits< double >::epsilon( ) );
            }

            // Check if expected state at end of integration step matches required state.
            // The reason this test has a different, higher tolerance than for a single integration
            // step is because the acummulated error builds up, since this test is based off of
            // Matlab time-series. The individual steps are accurate to the tolerance used in
            // Case 1.
            BOOST_CHECK_CLOSE_FRACTION( matlabDiscreteEventBenchmarkData( i, STATE_COLUMN_INDEX ),
                                        integrator.getCurrentState( )( 0 ),
                                        fullIntegrationTestTolerance );

            // Check if a discrete event is schedule to take place, and execute discrete event
            // if so.
            if ( std::fabs( integrator.getCurrentIndependentVariable( )
                            - TIME_OF_DISCRETE_EVENT_1 )
                 < std::numeric_limits< double >::epsilon( )
                 || std::fabs( integrator.getCurrentIndependentVariable( )
                               - TIME_OF_DISCRETE_EVENT_2 )
                    < std::numeric_limits< double >::epsilon( )
                 || std::fabs( integrator.getCurrentIndependentVariable( )
                               - TIME_OF_DISCRETE_EVENT_3 )
                    < std::numeric_limits< double >::epsilon( ) )
            {
                // Increment loop to next row in matrix, which contains the discrete event that
                // affects the state.
                i++;

                // Modify the current state based on the discrete event
                integrator.modifyCurrentState(
                            Eigen::VectorXd::Constant(
                                1, matlabDiscreteEventBenchmarkData( i, STATE_COLUMN_INDEX ) ) );

                // Check that it is now not possible to roll back.
                BOOST_CHECK( !integrator.rollbackToPreviousState( ) );
            }
        }

        // Check that it is possible to roll back to the previous step. This is allowable, because
        // the last action performed by the integrator is the performIntegrationStep() function.
        BOOST_CHECK( integrator.rollbackToPreviousState( ) );

        // Check that the rolled back time is as required. This test should be exact.
        BOOST_CHECK_EQUAL( lastTime, integrator.getCurrentIndependentVariable( ) );

        // Check that the rolled back state is as required. This test should be exact.
        BOOST_CHECK_EQUAL( lastState( 0 ), integrator.getCurrentState( )( 0 ) );

        // Check that it is now not possible to roll back to the previous step.
        BOOST_CHECK( !integrator.rollbackToPreviousState( ) );
    }
}

} // namespace matlab_numerical_integrator_tests
} // namespace unit_tests
} // namespace tudat

#endif // TUDAT_MATLAB_NUMERICAL_INTEGRATOR_TEST_H

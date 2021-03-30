/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *     If you are unable to replicate the error control mechanism for adaptive step size methods,
 *     you can force the integrator to accept all steps taken based on the time-data in the 
 *     benchmark data files. The easiest way to do this is to set the relative and absolute
 *     tolerances to a large value, so that in effect the error control mechanism doesn't kick
 *     in.
 *
 */

#ifndef TUDAT_NUMERICAL_INTEGRATOR_TESTS_H
#define TUDAT_NUMERICAL_INTEGRATOR_TESTS_H

#include <Eigen/Core>

#include "tudat/math/integrators/numericalIntegrator.h"
#include "tudat/math/integrators/reinitializableNumericalIntegrator.h"

namespace tudat
{
namespace unit_tests
{
namespace numerical_integrator_tests
{

//! Set column indices in matrices.
const int TIME_COLUMN_INDEX = 0;
const int STATE_COLUMN_INDEX = 1;
const int FIRST_ROW = 0;
const int SECOND_ROW = 1;

//! Use integrateTo() to integrate one step forward in time.
/*!
 * Uses the integrateTo() function provided with Tudat integrators to make one step forward in
 * time.
 * \param benchmarkData Benchmark integration data.
 * \param singleStepTestTolerance Tolerance used to compare computed final state with benchmark
 *          data.
 * \param integrator Shared-pointer to numerical integrator used.
 * \sa NumericalIntegrator.
 */
inline void executeOneIntegrateToStep(
        const Eigen::MatrixXd benchmarkData,
        const double singleStepTestTolerance,
        const numerical_integrators::NumericalIntegratorXdPointer integrator )
{
    // Use integrateTo() to integrate one step forward in time and check results against benchmark
    // data.

    // Check that it is now not possible to roll back.
    BOOST_CHECK( !integrator->rollbackToPreviousState( ) );

    // Use integrateTo() to integrate one step.
    integrator->integrateTo( benchmarkData( SECOND_ROW, TIME_COLUMN_INDEX ),
                             benchmarkData( SECOND_ROW, TIME_COLUMN_INDEX )
                             - benchmarkData( FIRST_ROW, TIME_COLUMN_INDEX ) );

    // Check if the expected time at the end of the integration step matches the required
    // time. This test should be exact.
    BOOST_CHECK_EQUAL( benchmarkData( SECOND_ROW, TIME_COLUMN_INDEX ),
                       integrator->getCurrentIndependentVariable( ) );

    // Check if expected state at end of integration step matches required state.
    BOOST_CHECK_CLOSE_FRACTION( benchmarkData( SECOND_ROW, STATE_COLUMN_INDEX ),
                                integrator->getCurrentState( )( 0 ), singleStepTestTolerance );

    // Roll back to the previous step. This should be possible since the integrateTo() function
    // was called above.
    BOOST_CHECK( integrator->rollbackToPreviousState( ) );

    // Check that the rolled back time is as required. This test should be exact.
    BOOST_CHECK_EQUAL( benchmarkData( FIRST_ROW, TIME_COLUMN_INDEX ),
                       integrator->getCurrentIndependentVariable( ) );

    // Check that the rolled back state is as required. This test should be exact.
    BOOST_CHECK_EQUAL( benchmarkData( FIRST_ROW, STATE_COLUMN_INDEX ),
                       integrator->getCurrentState( )( 0 ) );

    // Check that it is now not possible to roll back.
    BOOST_CHECK( !integrator->rollbackToPreviousState( ) );
}

//! Use performIntegrationStep() to integrate to specified time in multiple steps.
/*!
 * Uses the performIntegrationStep() function provided with Tudat integrations to integrate to a
 * specified time, with all of the intermediate steps checked against benchmark data.
 * \param benchmarkData Benchmark integration data.
 * \param singleStepTestTolerance Tolerance used to compare computed intermiedate states with
 *          benchmark data.
 * \param fullIntegrationTestTolerance Tolerance used to compare computed final state with
 *          benchmark data. This tolerance is greater than for a single step, as the errors
 *          accumulate during the integration.
 * \param integrator Shared-pointer to numerical integrator used.
 * \sa NumericalIntegrator.
 */
inline void performIntegrationStepToSpecifiedTime(
        const Eigen::MatrixXd benchmarkData,
        const double singleStepTestTolerance,
        const double fullIntegrationTestTolerance,
        const numerical_integrators::NumericalIntegratorXdPointer integrator )
{
    // Use performIntegrationstep() to integrate to specified time in multiple steps and check
    // results against benchmark data.

    // Declare last time and state.
    double lastTime = TUDAT_NAN;
    Eigen::VectorXd lastState = Eigen::VectorXd::Zero( 1 );

    // Loop through all the integration steps in the benchmark data.
    for ( int i = 1; i < benchmarkData.rows( ); i++ )
    {
        // Set the step size based on the benchmark data.
        double stepSize = benchmarkData( i, TIME_COLUMN_INDEX )
                - benchmarkData( i - 1, TIME_COLUMN_INDEX );

        // Store last time and state.
        lastTime = integrator->getCurrentIndependentVariable( );
        lastState = integrator->getCurrentState( );

        // Perform the integration step.
        integrator->performIntegrationStep( stepSize );

        // Check if the expected time at the end of the integration step matches the required
        // time. This test should be accurate to effectively machine precision.
        if ( std::fabs( integrator->getCurrentIndependentVariable( ) )
             < std::numeric_limits< double >::epsilon( ) )
        {
            BOOST_CHECK_SMALL( benchmarkData( i, TIME_COLUMN_INDEX ),
                               singleStepTestTolerance );
        }

        else
        {
            BOOST_CHECK_CLOSE_FRACTION( benchmarkData( i, TIME_COLUMN_INDEX ),
                                        integrator->getCurrentIndependentVariable( ),
                                        singleStepTestTolerance );
        }

        // Check if expected state at end of integration step matches required state.
        // The reason this test has a different, higher tolerance than for a single integration
        // step is because the acummulated error builds up, since this test is based off of
        // a time-series. The individual steps are accurate to singleStepTestTolerance.
        BOOST_CHECK_CLOSE_FRACTION( benchmarkData( i, STATE_COLUMN_INDEX ),
                                    integrator->getCurrentState( )( 0 ),
                                    fullIntegrationTestTolerance );
    }

    // Roll back to the previous step. This should be possible since the
    // performIntegrationStep() function was called above.
    BOOST_CHECK( integrator->rollbackToPreviousState( ) );

    // Check that the rolled back time is as required. This test should be exact.
    BOOST_CHECK_EQUAL( lastTime, integrator->getCurrentIndependentVariable( ) );

    // Check that the rolled back state is as required. This test should be exact.
    BOOST_CHECK_EQUAL( lastState( 0 ), integrator->getCurrentState( )( 0 ) );

    // Check that it is now not possible to roll back.
    BOOST_CHECK( !integrator->rollbackToPreviousState( ) );
}

//! Use integrateTo() to integrate to specified time in one step.
/*!
 * Uses the integrateTo() function provided with Tudat integrators to integrate to specified
 * time. All intermediate steps are used internally in the integrator, with the final step
 * being handled by the integrateTo() function.
 * \param benchmarkData Benchmark integration data.
 * \param fullIntegrationTestTolerance Tolerance used to compare computed final state with
 *          benchmark data. This tolerance is greater than for a single step, as the errors
 *          accumulate during the integration.
 * \param integrator Shared-pointer to numerical integrator used.
 * \param specifiedTime Time to integrator to.
 * \sa NumericalIntegrator.
 */
inline void executeIntegrateToToSpecifiedTime(
        const Eigen::MatrixXd benchmarkData,
        const double fullIntegrationTestTolerance,
        const numerical_integrators::NumericalIntegratorXdPointer integrator,
        const double specifiedTime )
{
    // Use integrateTo() to integrate to specified time and and check results against benchmark
    // data.

    // Set index of last row of matrix containing benchmark data.
    const int indexLastRowBenchmarkData = benchmarkData.rows( ) - 1;

    // Use integrateTo() to integrate to final time.
    integrator->integrateTo( specifiedTime,
                             benchmarkData( SECOND_ROW, TIME_COLUMN_INDEX )
                             - benchmarkData( FIRST_ROW, TIME_COLUMN_INDEX ) );

    // Check if the expected time at the end of the integration step matches the required
    // time. This test should be exact.
    BOOST_CHECK_EQUAL( benchmarkData( indexLastRowBenchmarkData, TIME_COLUMN_INDEX ),
                       integrator->getCurrentIndependentVariable( ) );

    // Check if expected state at end of integration step matches required state.
    // The reason this test has a different, higher tolerance than for a single integration
    // step is because the acummulated error builds up, since this test is based off of
    // a time-series. The individual steps are accurate to the tolerance used in
    // Case 1.
    BOOST_CHECK_CLOSE_FRACTION( benchmarkData( indexLastRowBenchmarkData,
                                               STATE_COLUMN_INDEX ),
                                integrator->getCurrentState( )( 0 ),
                                fullIntegrationTestTolerance );

    // Roll back to the previous step. This should be possible since the integrateTo() function
    // was called above.
    BOOST_CHECK( integrator->rollbackToPreviousState( ) );

    // Check that it is now not possible to roll back.
    BOOST_CHECK( !integrator->rollbackToPreviousState( ) );
}

//! Use performIntegrationStep() to integrate to specified time in multiple steps, including
//! discrete events.
/*!
 * Uses the performIntegrationStep() function provided with Tudat integrations to integrate to a
 * specified time, with all of the intermediate steps checked against benchmark data, and discrete
 * events, that instantly change the state, incorporated too.
 * \param benchmarkData Benchmark integration data.
 * \param singleStepTestTolerance Tolerance used to compare computed intermiedate states with
 *          benchmark data.
 * \param fullIntegrationTestTolerance Tolerance used to compare computed final state with
 *          benchmark data. This tolerance is greater than for a single step, as the errors
 *          accumulate during the integration.
 * \param integrator Shared-pointer to re-initializable numerical integrator used.
 * \sa ReinitializableNumericalIntegrator.
 */
inline void performIntegrationStepToSpecifiedTimeWithEvents(
        const Eigen::MatrixXd benchmarkData,
        const double singleStepTestTolerance,
        const double fullIntegrationTestTolerance,
        const numerical_integrators::ReinitializableNumericalIntegratorXdPointer integrator )
{
    // Use performIntegrationstep() to integrate to specified time in multiple steps, including
    // discrete events and check results against benchmark data. The discrete events are enabled
    // by modifying the state instantaneously. These discrete events cannot be rolled back.

    // Declare last time and state.
    double lastTime = TUDAT_NAN;
    Eigen::VectorXd lastState = Eigen::VectorXd::Zero( 1 );

    // Loop through all the integration steps in the benchmark data.
    for ( int i = 1; i < benchmarkData.rows( ); i++ )
    {
        // Set the step size based on the benchmark data generated.
        const double stepSize = benchmarkData( i, TIME_COLUMN_INDEX )
                - benchmarkData( i - 1, TIME_COLUMN_INDEX );

        // Store last time and state.
        lastTime = integrator->getCurrentIndependentVariable( );
        lastState = integrator->getCurrentState( );

        // Perform the integration step.
        integrator->performIntegrationStep( stepSize );

        // Check if the expected time at the end of the integration step matches the required
        // time. This test should be accurate to effectively machine precision.
        if ( std::fabs( benchmarkData( i, TIME_COLUMN_INDEX ) )
             < std::numeric_limits< double >::min( ) )
        {
            BOOST_CHECK_SMALL( integrator->getCurrentIndependentVariable( ),
                               std::numeric_limits< double >::epsilon( ) );
        }

        else
        {
            BOOST_CHECK_CLOSE_FRACTION(
                        benchmarkData( i, TIME_COLUMN_INDEX ),
                        integrator->getCurrentIndependentVariable( ),
                        singleStepTestTolerance );
        }

        // Check if expected state at end of integration step matches required state.
        // The reason this test has a different, higher tolerance than for a single integration
        // step is because the acummulated error builds up, since this test is based off of
        // a time-series. The individual steps are accurate to the tolerance used in
        // Case 1.
        BOOST_CHECK_CLOSE_FRACTION( benchmarkData( i, STATE_COLUMN_INDEX ),
                                    integrator->getCurrentState( )( 0 ),
                                    fullIntegrationTestTolerance );

        // Check if a discrete event is schedule to take place, and execute discrete event
        // if so.
        if ( i < benchmarkData.rows( ) - 1
             && std::fabs( benchmarkData( i + 1, TIME_COLUMN_INDEX )
                           - benchmarkData( i, TIME_COLUMN_INDEX ) )
             < std::numeric_limits< double >::epsilon( ) )
        {
            // Increment loop to next row in matrix, which contains the discrete event that
            // affects the state.
            i++;

            // Modify the current state based on the discrete event
            integrator->modifyCurrentState(
                        Eigen::VectorXd::Constant(
                            1, benchmarkData( i, STATE_COLUMN_INDEX ) ) );

            // Check that it is now not possible to roll back.
            BOOST_CHECK( !integrator->rollbackToPreviousState( ) );
        }
    }

    // Check that it is possible to roll back to the previous step. This is allowable, because
    // the last action performed by the integrator is the performIntegrationStep() function.
    BOOST_CHECK( integrator->rollbackToPreviousState( ) );

    // Check that the rolled back time is as required. This test should be exact.
    BOOST_CHECK_EQUAL( lastTime, integrator->getCurrentIndependentVariable( ) );

    // Check that the rolled back state is as required. This test should be exact.
    BOOST_CHECK_EQUAL( lastState( 0 ), integrator->getCurrentState( )( 0 ) );

    // Check that it is now not possible to roll back to the previous step.
    BOOST_CHECK( !integrator->rollbackToPreviousState( ) );
}
} // namespace numerical_integrator_tests
} // namespace unit_tests
} // namespace tudat

#endif // TUDAT_NUMERICAL_INTEGRATOR_TESTS_H

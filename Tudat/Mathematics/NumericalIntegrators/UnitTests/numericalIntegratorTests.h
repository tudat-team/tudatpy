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
 *      120329    K. Kumar          File created.
 *      120404    K. Kumar          Added unit test for discrete events.
 *      121106    K. Kumar          Generalized code to work with any numerical integrator.
 *      121218    K. Kumar          Migrated to Tudat Core.
 *      130115    K. Kumar          Finalized tests, merged forwards and backwards integration
 *                                  test functions.
 *
 *    References
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

#include "Tudat/Mathematics/NumericalIntegrators/numericalIntegrator.h"
#include "Tudat/Mathematics/NumericalIntegrators/reinitializableNumericalIntegrator.h"

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
void executeOneIntegrateToStep(
        const Eigen::MatrixXd benchmarkData,
        const double singleStepTestTolerance,
        const numerical_integrators::NumericalIntegratorXdPointer integrator );

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
void performIntegrationStepToSpecifiedTime(
        const Eigen::MatrixXd benchmarkData,
        const double singleStepTestTolerance,
        const double fullIntegrationTestTolerance,
        const numerical_integrators::NumericalIntegratorXdPointer integrator );

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
void executeIntegrateToToSpecifiedTime(
        const Eigen::MatrixXd benchmarkData,
        const double fullIntegrationTestTolerance,
        const numerical_integrators::NumericalIntegratorXdPointer integrator,
        const double specifiedTime );

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
void performIntegrationStepToSpecifiedTimeWithEvents(
        const Eigen::MatrixXd benchmarkData,
        const double singleStepTestTolerance,
        const double fullIntegrationTestTolerance,
        const numerical_integrators::ReinitializableNumericalIntegratorXdPointer integrator );

} // namespace numerical_integrator_tests
} // namespace unit_tests
} // namespace tudat

#endif // TUDAT_NUMERICAL_INTEGRATOR_TESTS_H

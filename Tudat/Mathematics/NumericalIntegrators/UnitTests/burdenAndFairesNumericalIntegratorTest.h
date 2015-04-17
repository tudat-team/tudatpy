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
 *      120404    K. Kumar          Added missing comments.
 *
 *    References
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole, 2001.
 *
 *    Notes
 *      Currently, this file is only used to test the RKF45 integrator. There is however more test
 *      data available in Burden and Faires (2001), for other integrators, that can be implemented
 *      for unit testing.
 *
 */

#ifndef TUDAT_BURDEN_AND_FAIRES_NUMERICAL_INTEGRATOR_TEST_H
#define TUDAT_BURDEN_AND_FAIRES_NUMERICAL_INTEGRATOR_TEST_H

#include <cmath>
#include <utility>

#include <Eigen/Core>

#include "Tudat/Basics/utilityMacros.h"
#include "Tudat/Mathematics/NumericalIntegrators/UnitTests/numericalIntegratorTestFunctions.h"

namespace tudat
{
namespace unit_tests
{

//! Class for Burden and Faires numerical integrators tests.
/*!
 * Class with function to compute new step size for Runge-Kutta variable step size integrators.
 * The function is based on the error estimation method given in Algorithm 5.3 in
 * (Burden and Faires, 2001).
 */
class BurdenAndFairesNumericalIntegratorTest
{
public:

    //! Default constructor.
    /*!
     * Default constructor, initializing input parameters.
     */
    BurdenAndFairesNumericalIntegratorTest( ) :
        relativeError_( Eigen::VectorXd::Zero( 1 ) ),
        lowerOrderEstimate_( Eigen::VectorXd::Zero( 1 ) ),
        higherOrderEstimate_( Eigen::VectorXd::Zero( 1 ) )
    { }

    //! Compute new step size using method taken from (Burden and Faires, 2001).
    /*!
     * Computes the new step size based on the error estimation method given in Algorithm 5.3 in
     * (Burden and Faires, 2001).
     */
    std::pair< double, bool > computeNewStepSize( const double stepSize, const double lowerOrder,
                                                  const double higherOrder,
                                                  const double safetyFactorForNextStepSize,
                                                  const Eigen::VectorXd& relativeErrorTolerance,
                                                  const Eigen::VectorXd& absoluteErrorTolerance,
                                                  const Eigen::VectorXd& lowerOrderEstimate,
                                                  const Eigen::VectorXd& higherOrderEstimate )
    {
        TUDAT_UNUSED_PARAMETER( higherOrder );
        TUDAT_UNUSED_PARAMETER( relativeErrorTolerance );

        // Set higher and lower order estimates.
        higherOrderEstimate_ = higherOrderEstimate;
        lowerOrderEstimate_ = lowerOrderEstimate;

        // Compute truncation error estimate based on higher and lower estimates.
        const Eigen::VectorXd trucationError = ( higherOrderEstimate
                                                 - lowerOrderEstimate ).array( ).abs( );

        // Compute relative error.
        relativeError_ = trucationError.array( ).abs( ) / stepSize;

        // Compute new step size for next integration step. (In case the step size for this
        // current integration step doesn't satisfy the tolerances set, this new step size is
        // used to recompute the current step.
        const double newStepSize = safetyFactorForNextStepSize * stepSize * std::pow(
                    absoluteErrorTolerance.array( ).minCoeff( )
                    / relativeError_.array( ).abs( ).maxCoeff( ), 1.0 / lowerOrder );

        // Check if the current integration step is accepted, based on the allowed absolute
        // tolerance.
        const bool isIntegrationStepAccepted = relativeError_.array( ).abs( ).maxCoeff( )
                < absoluteErrorTolerance.array( ).minCoeff( );

        // Return the computed new step size and a flag whether the current step is accepted or
        // not.
        return std::make_pair( newStepSize, isIntegrationStepAccepted );
    }

    //! Compute state derivative for Burden and Faires example.
    /*!
     * Computes state derivative. The state derivative function defined corresponds to
     * Example 3, pg. 278 in (Burden and Faires, 2001).
     * The initial-value problem is:
     *
     * \f[
     *      y' = y - t^{ 2 } + 1
     * \f]
     *
     * with \f$ 0 \leq t \leq 2 \f$ and \f$ y( 0 ) = 0.5 \f$.
     * \param time Time at which the state derivative needs to be evaluated.
     * \param state State, length of which should be 1, at which the state derivative needs to
     *          be evaluated.
     * \return State derivative, length equal to state values according to above expression.
     */
    Eigen::VectorXd computeStateDerivative( const double time, const Eigen::VectorXd& state )
    {
        return numerical_integrator_test_functions::
                computeNonAutonomousModelStateDerivative( time, state );
    }

    //! Relative error.
    /*!
     * Relative error computed by computeNewStepSize( ).
     */
    Eigen::VectorXd relativeError_;

    //! Lower order estimate.
    /*!
     * Lower order estimate computed for error control.
     */
    Eigen::VectorXd lowerOrderEstimate_;

    //! Higher order estimate.
    /*!
     * Higher order estimate computed for error control.
     */
    Eigen::VectorXd higherOrderEstimate_;

protected:

private:
};

} // namespace unit_tests
} // namespace tudat

#endif // TUDAT_BURDEN_AND_FAIRES_NUMERICAL_INTEGRATOR_TEST_H

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
 *      120329    K. Kumar          File created.
 *      120404    K. Kumar          Added missing comments.
 *
 *    References
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole, 2001.
 *
 */

#ifndef TUDAT_BURDEN_AND_FAIRES_NUMERICAL_INTEGRATOR_TEST_H
#define TUDAT_BURDEN_AND_FAIRES_NUMERICAL_INTEGRATOR_TEST_H

#include <cmath>
#include <utility>

#include <Eigen/Core>

#include <TudatCore/Basics/utilityMacros.h>

#include <Tudat/Mathematics/NumericalIntegrators/UnitTests/numericalIntegratorTestFunctionSuite.h>

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
        return numerical_integrator_test_function_suite::computeNonAutonomousModelStateDerivative(
                    time, state );
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

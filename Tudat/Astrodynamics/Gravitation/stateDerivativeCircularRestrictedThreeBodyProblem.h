/*    Copyright (c) 2010 Delft University of Technology.
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
 *      110905    L. van der Ham    First creation of code.
 *      110922    L. van der Ham    Added computation of Jacobian energy constant.
 *      111107    L. van der Ham    Removed option planar.
 *      111110    K. Kumar          Minor comment and naming modifications.
 *      120221    L. van der Ham    Removed computation of Jacobian energy constant.
 *      120309    K. Kumar          Updated code to latest Tudat standards; replaced set-function
 *                                  for mass parameter with constructor input; renamed class and
 *                                  file.
 *
 *    References
 *        Wakker, K.F., "Astrodynamics I, AE4-874", Delft University of Technology, 2007.
 *
 */

#ifndef TUDAT_STATE_DERIVATIVE_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H
#define TUDAT_STATE_DERIVATIVE_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H

#include <Eigen/Core>

namespace tudat
{
namespace astrodynamics
{
namespace gravitation
{
namespace circular_restricted_three_body_problem
{

//! State elements indices in CRTBP.
/*!
 * State elements indices in CRTBP. These are indices for Cartesian elements in normalized units,
 * in a co-rotating reference frame (co-rotating with the primaries).
 */
enum StateElementIndices
{
    normalizedXPositionIndex,
    normalizedYPositionIndex,
    normalizedZPositionIndex,
    normalizedXVelocityIndex,
    normalizedYVelocityIndex,
    normalizedZVelocityIndex
};

class StateDerivativeCircularRestrictedThreeBodyProblem
{
public:

    //! Default constructor.
    /*!
     * Default constructor that defines the state derivative for a given CRTBP system.
     * \param massParameter Mass parameter of CRTBP.
     */
    StateDerivativeCircularRestrictedThreeBodyProblem( double massParameter )
        : massParameter_ ( massParameter )
    { }

    //! Compute state derivative.
    /*!
     * Computes the state derivative of CRTBP.
     * \param time Time.
     * \param cartesianState Cartesian state.
     * \return State derivative.
     */
    Eigen::VectorXd computeStateDerivative( const double time,
                                            const Eigen::VectorXd& cartesianState );

protected:

private:

    //! Mass parameter.
    /*!
     * Value of mass parameter for the CRTBP.
     */
    double massParameter_;
};

} // namespace circular_restricted_three_body_problem
} // namespace gravitation
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_STATE_DERIVATIVE_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H

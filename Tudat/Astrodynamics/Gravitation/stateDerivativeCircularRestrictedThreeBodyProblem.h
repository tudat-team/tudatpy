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
 *      110905    L. van der Ham    File created.
 *      110922    L. van der Ham    Added computation of Jacobian energy constant.
 *      111107    L. van der Ham    Removed option planar.
 *      111110    K. Kumar          Minor comment and naming modifications.
 *      120221    L. van der Ham    Removed computation of Jacobian energy constant.
 *      120309    K. Kumar          Updated code to latest Tudat standards; replaced set-function
 *                                  for mass parameter with constructor input; renamed class and
 *                                  file.
 *      120426    K. Kumar          Added enum for state derivative acceleration elements.
 *      130121    K. Kumar          Added shared-ptr typedef; updated VectorXd to Vector6d.
 *
 *    References
 *        Wakker, K.F., "Astrodynamics I, AE4-874", Delft University of Technology, 2007.
 *
 *    Notes
 *
 */

#ifndef TUDAT_STATE_DERIVATIVE_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H
#define TUDAT_STATE_DERIVATIVE_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

namespace tudat
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
    xCartesianPositionIndex,
    yCartesianPositionIndex,
    zCartesianPositionIndex,
    xCartesianVelocityIndex,
    yCartesianVelocityIndex,
    zCartesianVelocityIndex
};

//! State derivative elements indices in CRTBP.
/*!
 * State derivative elements indices in CRTBP. These are indices for Cartesian elements in
 * normalized units, in a co-rotating reference frame (co-rotating with the primaries).
 */
enum StateDerivativeElementIndices
{
    xAccelerationIndex = 3,
    yAccelerationIndex = 4,
    zAccelerationIndex = 5
};

//! State derivative model class for CRTBP.
/*!
 * Class that contains the state derivative model for the CRTBP.
 */
class StateDerivativeCircularRestrictedThreeBodyProblem
{
public:

    //! Default constructor.
    /*!
     * Default constructor that defines the state derivative for a given CRTBP system.
     * \param aMassParameter A value for mass parameter of CRTBP.
     */
    StateDerivativeCircularRestrictedThreeBodyProblem( const double aMassParameter )
        : massParameter ( aMassParameter )
    { }

    //! Compute state derivative.
    /*!
     * Computes the state derivative of CRTBP.
     * \param time Time.
     * \param cartesianState Cartesian state.
     * \return State derivative.
     */
    Eigen::Vector6d computeStateDerivative(
            const double time, const Eigen::Vector6d& cartesianState );

protected:

private:

    //! Mass parameter.
    /*!
     * Value of mass parameter for the CRTBP.
     */
    double massParameter;
};

//! Typedef for shared-pointer to StateDerivativeCircularRestrictedThreeBodyProblem object.
typedef boost::shared_ptr< StateDerivativeCircularRestrictedThreeBodyProblem >
StateDerivativeCircularRestrictedThreeBodyProblemPointer;

} // namespace circular_restricted_three_body_problem
} // namespace gravitation
} // namespace tudat

#endif // TUDAT_STATE_DERIVATIVE_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H

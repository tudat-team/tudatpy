/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *        Wakker, K.F., "Astrodynamics I, AE4-874", Delft University of Technology, 2007.
 *
 */

#ifndef TUDAT_STATE_DERIVATIVE_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H
#define TUDAT_STATE_DERIVATIVE_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H

#include <memory>

#include <Eigen/Core>

#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{
namespace propagators
{
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
typedef std::shared_ptr< StateDerivativeCircularRestrictedThreeBodyProblem >
StateDerivativeCircularRestrictedThreeBodyProblemPointer;

} // namespace propagators

} // namespace tudat

#endif // TUDAT_STATE_DERIVATIVE_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H

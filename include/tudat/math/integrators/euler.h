/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *
 */

#ifndef TUDAT_EULER_INTEGRATOR_H
#define TUDAT_EULER_INTEGRATOR_H

#include <memory>

#include <Eigen/Core>

#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/math/integrators/rungeKuttaFixedStepSizeIntegrator.h"

namespace tudat
{
namespace numerical_integrators
{

//! Class that implements the Euler integrator.
/*!
 * Class that implements the Euler, fixed order, fixed step size integrator.
 * \tparam StateType The type of the state. This type should support addition with
 *          StateDerivativeType
 * \tparam StateDerivativeType The type of the state derivative. This type should support
 *          multiplication with IndependentVariableType and doubles.
 * \tparam IndependentVariableType The type of the independent variable.
 * \sa NumericalIntegrator.
 */
template< typename IndependentVariableType = double, typename StateType = Eigen::VectorXd,
           typename StateDerivativeType = Eigen::VectorXd, typename TimeStepType = IndependentVariableType >
class EulerIntegrator
        : public numerical_integrators::RungeKuttaFixedStepSizeIntegrator<
        IndependentVariableType, StateType, StateDerivativeType, TimeStepType >
{
public:

    //! Typedef for the base class.
    /*!
     * Typedef of the base class with all template parameters filled in.
     */
    typedef numerical_integrators::RungeKuttaFixedStepSizeIntegrator<
    IndependentVariableType, StateType, StateDerivativeType, TimeStepType > RungeKuttaFixedStepSizeIntegratorBase;

    //! Typedef for the state derivative function.
    /*!
     * Typedef to the state derivative function inherited from the base class.
     */
    typedef typename RungeKuttaFixedStepSizeIntegratorBase::
    StateDerivativeFunction StateDerivativeFunction;

    //! Default constructor.
    /*!
     * Default constructor, taking a state derivative function as argument.
     * \param stateDerivativeFunction State derivative function.
     * \param intervalStart The start of the integration interval.
     * \param initialState The initial state.
     * \sa NumericalIntegrator::NumericalIntegrator.
     */
    EulerIntegrator( const StateDerivativeFunction& stateDerivativeFunction,
                     const IndependentVariableType intervalStart,
                     const StateType& initialState )
        : RungeKuttaFixedStepSizeIntegratorBase( stateDerivativeFunction, intervalStart, initialState, forwardEuler )
    {
    }

protected:

};

//! Typedef of Euler integrator (state/state derivative = VectorXd, independent variable = double).
/*!
 * Typedef of an Euler integrator with VectorXds as state and state derivative and double as
 * independent variable.
 */
typedef EulerIntegrator< > EulerIntegratorXd;

//! Typedef of a scalar Euler integrator.
/*!
 * Typedef of a Euler integrator with doubles as state and state derivative and independent variable.
 */
typedef EulerIntegrator< double, double, double > EulerIntegratord;

//! Typedef for a shared-pointer to default Euler integrator.
/*!
 * Typedef for a shared-pointer to an Euler integrator with VectorXds as state and state derivative and double
 * as independent variable.
 */
typedef std::shared_ptr< EulerIntegratorXd > EulerIntegratorXdPointer;

//! Typedef of pointer to a scalar Euler integrator.
/*!
 * Typedef of pointer to an Euler integrator with doubles as state and state derivative and
 * independent variable.
 */
typedef std::shared_ptr< EulerIntegratord > EulerIntegratordPointer;

} // namespace numerical_integrators
} // namespace tudat

#endif // TUDAT_EULER_INTEGRATOR_H

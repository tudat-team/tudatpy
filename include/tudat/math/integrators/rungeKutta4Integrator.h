/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */
#ifdef NDEBUG
#ifdef TUDAT_BUILD_GNU
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif // TUDAT_BUILD_GNU
#endif // NDEBUG

#ifndef TUDAT_RUNGE_KUTTA_4_INTEGRATOR_H
#define TUDAT_RUNGE_KUTTA_4_INTEGRATOR_H

#include <memory>

#include <Eigen/Core>

#include "tudat/math/integrators/reinitializableNumericalIntegrator.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/math/integrators/rungeKuttaFixedStepSizeIntegrator.h"

namespace tudat
{

namespace numerical_integrators
{

//! Class that implements the Runge-Kutta 4 integrator.
/*!
 * Class that implements the Runge-Kutta 4, fixed order, fixed step size integrator.
 * \tparam StateType The type of the state. This type should support addition with
 *          StateDerivativeType.
 * \tparam StateDerivativeType The type of the state derivative. This type should support
 *          multiplication with IndependentVariableType and doubles.
 * \tparam IndependentVariableType The type of the independent variable.
 * \sa NumericalIntegrator.
 */
template< typename IndependentVariableType = double,
         typename StateType = Eigen::VectorXd,
           typename StateDerivativeType = Eigen::VectorXd,
         typename TimeStepType = IndependentVariableType >
class RungeKutta4Integrator
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
     */
    RungeKutta4Integrator( const StateDerivativeFunction& stateDerivativeFunction,
                           const IndependentVariableType intervalStart,
                           const StateType& initialState,
                           const TimeStepType& stepSize )
        : RungeKuttaFixedStepSizeIntegratorBase( stateDerivativeFunction, intervalStart, initialState, stepSize, rungeKutta4Classic )
    {
    }

protected:

};

//extern template class RungeKutta4Integrator < double, Eigen::VectorXd, Eigen::VectorXd >;
//extern template class RungeKutta4Integrator < double, Eigen::Vector6d, Eigen::Vector6d >;
//extern template class RungeKutta4Integrator < double, Eigen::MatrixXd, Eigen::MatrixXd >;


//! Typedef of RK4 integrator (state/state derivative = VectorXd, independent variable = double).
/*!
 * Typedef of a RK4 integrator with VectorXds as state and state derivative and double as
 * independent variable.
 */
typedef RungeKutta4Integrator< > RungeKutta4IntegratorXd;

//! Typedef of a scalar RK4 integrator.
/*!
 * Typedef of an RK4 integrator with doubles as state and state derivative and independent
 * variable.
 */
typedef RungeKutta4Integrator< double, double, double > RungeKutta4Integratord;

//! Typedef of pointer to default RK4 integrator
/*!
 * Typedef of pointer to a RK4 integrator with VectorXds as state and state derivative and double
 * as independent variable.
 */
typedef std::shared_ptr< RungeKutta4IntegratorXd > RungeKutta4IntegratorXdPointer;

//! Typedef of pointer to a scalar RK4 integrator.
/*!
 * Typedef of pointer to an RK4 integrator with doubles as state and state derivative and
 * independent variable.
 */
typedef std::shared_ptr< RungeKutta4Integratord > RungeKutta4IntegratordPointer;

} // namespace numerical_integrators

} // namespace tudat

#endif // TUDAT_RUNGE_KUTTA_4_INTEGRATOR_H

#ifdef NDEBUG
#ifdef TUDAT_BUILD_GNU
// turn the warnings back on
#pragma GCC diagnostic pop
#endif // TUDAT_BUILD_GNU
#endif // NDEBUG

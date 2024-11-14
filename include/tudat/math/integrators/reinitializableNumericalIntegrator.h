/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_REINITIALIZABLE_NUMERICAL_INTEGRATOR_H
#define TUDAT_REINITIALIZABLE_NUMERICAL_INTEGRATOR_H

#include <Eigen/Core>

#include "tudat/math/integrators/numericalIntegrator.h"

namespace tudat
{

namespace numerical_integrators
{

//! Base class for reinitializable numerical integrators.
/*!
 * This base class enables derived classes to reinitialize the current state of the numerical
 * integrator.
 * \tparam IndependentVariableType The type of the independent variable.
 * \tparam StateType The type of the state. This type should support addition with
 *          StateDerivativeType.
 * \tparam StateDerivativeType The type of the state derivative. This type should support
 *          multiplication with IndependentVariableType and doubles.
 * \sa NumericalIntegrator.
 */
template< typename IndependentVariableType = double, typename StateType = Eigen::VectorXd,
           typename StateDerivativeType = Eigen::VectorXd, typename TimeStepType = IndependentVariableType >
class ReinitializableNumericalIntegrator :
        public numerical_integrators::NumericalIntegrator<
        IndependentVariableType, StateType, StateDerivativeType, TimeStepType >
{
protected:

    //! Typedef of the base class.
    /*!
     * Typedef of the base class with all template parameters filled in.
     */
    typedef numerical_integrators::NumericalIntegrator<
    IndependentVariableType, StateType, StateDerivativeType, TimeStepType > NumericalIntegratorBase;

public:

    //! Default constructor.
    ReinitializableNumericalIntegrator(
            const typename NumericalIntegratorBase::StateDerivativeFunction& aStateDerivativeFunction )
        : NumericalIntegratorBase( aStateDerivativeFunction )
    { }

    //! Default destructor.
    virtual ~ReinitializableNumericalIntegrator( ) { }

protected:

private:

};

//extern template class ReinitializableNumericalIntegrator < double, Eigen::VectorXd, Eigen::VectorXd >;
//extern template class ReinitializableNumericalIntegrator < double, Eigen::Vector6d, Eigen::Vector6d >;
//extern template class ReinitializableNumericalIntegrator < double, Eigen::MatrixXd, Eigen::MatrixXd >;

//! Typedef for shared-pointer to default, re-initializable numerical integrator.
/*!
 * Typedef for shared-pointer to a default, re-initializable numerical integrator
 * (IndependentVariableType = double, StateType = Eigen::VectorXd,
 * StateDerivativeType = Eigen::VectorXd).
 */
typedef std::shared_ptr< ReinitializableNumericalIntegrator< > >
ReinitializableNumericalIntegratorXdPointer;

//! Typedef for a shared-pointer to a scalar, re-initializable numerical integrator.
/*!
 * Typedef for shared-pointer to a scalar numerical, re-initializable integrator
 * (IndependentVariableType = double, StateType = double, StateDerivativeType = double).
 */
typedef std::shared_ptr< ReinitializableNumericalIntegrator< double, double, double > >
ReinitializableNumericalIntegratordPointer;

} // namespace numerical_integrators

} // namespace tudat

#endif // TUDAT_REINITIALIZABLE_NUMERICAL_INTEGRATOR_H

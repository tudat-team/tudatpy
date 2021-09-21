/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROPAGATECOVARIANCE_H
#define TUDAT_PROPAGATECOVARIANCE_H

#include<tudat/astro/propagators/stateTransitionMatrixInterface.h>

namespace tudat
{

namespace propagators
{

//! Function to retrueve full state transition and sensitivity matrices at epochs
/*!
 *  Function to retrueve full state transition and sensitivity matrices at epochs from associated interface object
 * \param fullVariationalEquationsSolutionHistory List of combined state transition and sensitivity matrices, given
 * at the epochs where the state covariance is to be computed
 * \param stateTransitionInterface Object that is used to obtain state transition and sensitivity matrices
 * \param evaluationTimes Times at which the variational equations are to be evaluated
 */
void getFullVariationalEquationsSolutionHistory(
        std::map< double, Eigen::MatrixXd >& fullVariationalEquationsSolutionHistory,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const std::vector< double > evaluationTimes );


//! Function to retrueve full state transition and sensitivity matrices at epochs
/*!
 *  Function to retrueve full state transition and sensitivity matrices at epochs from associated interface object
 * \param fullVariationalEquationsSolutionHistory List of combined state transition and sensitivity matrices, given
 * at the epochs where the state covariance is to be computed
 * \param stateTransitionInterface Object that is used to obtain state transition and sensitivity matrices
 * \param timeStep Time step with which full state transition and sensitivity matrices is to be provided
 * \param initialTime Initial time at which full state transition and sensitivity matrices is to be provided
 * \param finalTime Final time at which full state transition and sensitivity matrices is to be provided
 */
void getFullVariationalEquationsSolutionHistory(
        std::map< double, Eigen::MatrixXd >& fullVariationalEquationsSolutionHistory,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const double timeStep,
        const double initialTime,
        const double finalTime );

//! Function to propagate full covariance at the initial time to state covariance at later times
/*!
 * Function to propagate full covariance at the initial time to state covariance at later times
 * \param propagatedCovariance List of state covariances at epochs (returned by reference)
 * \param initialCovariance Full covariance at initial time
 * \param fullVariationalEquationsSolutionHistory List of combined state transition and sensitivity matrices, given
 * at the epochs where the state covariance is to be computed
 */
void propagateCovariance(
        std::map< double, Eigen::MatrixXd >& propagatedCovariance,
        const Eigen::MatrixXd& initialCovariance,
        const std::map< double, Eigen::MatrixXd >& fullVariationalEquationsSolutionHistory );

//! Function to propagate full covariance at the initial time to state covariance at later times
/*!
 * Function to propagate full covariance at the initial time to state covariance at later times
 * \param propagatedCovariance List of state covariances at epochs (returned by reference)
 * \param initialCovariance Full covariance at initial time
 * \param stateTransitionInterface Object that is used to obtain state transition and sensitivity matrices
 * \param evaluationTimes Times at which the covariance is to be evaluated
 */
void propagateCovariance(
        std::map< double, Eigen::MatrixXd >& propagatedCovariance,
        const Eigen::MatrixXd& initialCovariance,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const std::vector< double > evaluationTimes );

//! Function to propagate full covariance at the initial time to state formal errors at later times
std::map< double, Eigen::MatrixXd > propagateCovariance(
        const Eigen::MatrixXd initialCovariance,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const std::vector< double > evaluationTimes );


//! Function to propagate full covariance at the initial time to state covariance at later times
/*!
 * Function to propagate full covariance at the initial time to state covariance at later times
 * \param propagatedCovariance List of state covariances at epochs (returned by reference)
 * \param initialCovariance Full covariance at initial time
 * \param stateTransitionInterface Object that is used to obtain state transition and sensitivity matrices
 * \param timeStep Time step with which covariance is to be provided
 * \param initialTime Initial time at which covariance is to be provided
 * \param finalTime Final time at which covariance is to be provided
 */
void propagateCovariance(
        std::map< double, Eigen::MatrixXd >& propagatedCovariance,
        const Eigen::MatrixXd& initialCovariance,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const double timeStep,
        const double initialTime,
        const double finalTime );

//! Function to propagate full covariance at the initial time to state formal errors at later times
/*!
 * Function to propagate full covariance at the initial time to state formal errors at later times
 * \param propagatedFormalErrors List of state formal errors at epochs (returned by reference)
 * \param initialCovariance Full covariance at initial time
 * \param stateTransitionInterface Object that is used to obtain state transition and sensitivity matrices
 * \param evaluationTimes Times at which the covariance is to be evaluated
 */
void propagateFormalErrors(
        std::map< double, Eigen::VectorXd >& propagatedFormalErrors,
        const Eigen::MatrixXd& initialCovariance,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const std::vector< double > evaluationTimes );

//! Function to propagate full covariance at the initial time to state formal errors at later times
std::map< double, Eigen::VectorXd > propagateFormalErrors(
        const Eigen::MatrixXd initialCovariance,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const std::vector< double > evaluationTimes );

//! Function to propagate full covariance at the initial time to state formal errors at later times
/*!
 * Function to propagate full covariance at the initial time to state formal errors at later times
 * \param propagatedFormalErrors List of state formal errors at epochs (returned by reference)
 * \param initialCovariance Full covariance at initial time
 * \param stateTransitionInterface Object that is used to obtain state transition and sensitivity matrices
 * \param timeStep Time step with which formal errors is to be provided
 * \param initialTime Initial time at which formal errors is to be provided
 * \param finalTime Final time at which formal errors is to be provided
 */
void propagateFormalErrors(
        std::map< double, Eigen::VectorXd >& propagatedFormalErrors,
        const Eigen::MatrixXd& initialCovariance,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const double timeStep,
        const double initialTime,
        const double finalTime );

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATECOVARIANCE_H

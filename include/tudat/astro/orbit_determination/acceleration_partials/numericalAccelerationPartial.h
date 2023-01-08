/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NUMERICALACCELERATIONPARTIAL_H
#define TUDAT_NUMERICALACCELERATIONPARTIAL_H

#include <functional>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/massRateModel.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/basic_astro/torqueModel.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/basics/basicTypedefs.h"

namespace tudat
{

namespace acceleration_partials
{

//! Dummy function used for update, performs no calculations.
/*!
 *  Dummy function used for update, performs no calculations.
 */
void emptyFunction( );

//! Dummy function used for update, performs no calculations.
/*!
 *  Dummy function used for update, performs no calculations.
 *  \param time Input parameter, not used by function.
 */
void emptyTimeFunction( const double time );

//! Function to numerical compute the partial derivative of an acceleration w.r.t. a body state.
/*!
 * Function to numerical compute the partial derivative of an acceleration w.r.t. a body state (position or velocity),
 * using a first-order central difference method.
 * \param setBodyState Function to reset the current state w.r.t. which the partial is to be computed.
 * \param accelerationModel Acceleration model for which the partial derivative is to be computed.
 * \param originalState Nominal state at which the partial derivative is to be computed.
 * \param statePerturbation Perturbation to the position or velocity that is to be used
 * \param startIndex Start index in the state vector from where the statePerturbation is to be added (i.e. 0 if
 * function is to compute partial w.r.t. position, 3 if w.r.t velocity).
 * \param updateFunction Function to update the required environment models following the change of the body state.
 * \param evaluationTime Time at which partial is to be evaluated (default NaN).
 * \return Numerical partial of the acceleration w.r.t. position or velocity (depending on function input).
 */
Eigen::Matrix3d calculateAccelerationWrtStatePartials(
        std::function< void( Eigen::Vector6d ) > setBodyState,
        std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        Eigen::Vector6d originalState,
        Eigen::Vector3d statePerturbation,
        int startIndex,
        std::function< void( ) > updateFunction = emptyFunction,
        const double evaluationTime = TUDAT_NAN );

Eigen::Vector3d calculateAccelerationWrtMassPartials(
        std::function< void( double ) > setBodyMass,
        std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        double originalMass,
        double massPerturbation,
        std::function< void( ) > updateFunction = emptyFunction,
        const double evaluationTime = TUDAT_NAN );


//! Function to numerical compute the partial derivative of a torque w.r.t. a body rotational state.
/*!
 * Function to numerical compute the partial derivative of an torque w.r.t. a body rotational state elements
 * using a first-order central difference method.
 * \param setBodyRotationalState Function to reset the current rotational state w.r.t. which the partial is to be computed.
 * \param torqueModel Torque model for which the partial derivative is to be computed.
 * \param originalRotationalState Nominal rotational state at which the partial derivative is to be computed.
 * \param statePerturbations Perturbation to the states that is to be used
 * \param startIndex Start index in the state vector from where the statePerturbation is to be added (i.e. 0 if
 * function is to compute partial w.r.t. quaternion, 4 if w.r.t angular velocity).
 * \param numberOfEntries Integer denoting the number of entries.
 * \param updateFunction Function to update the required environment models following the change of the body state.
 * \param evaluationTime Time at which partial is to be evaluated (default NaN).
 * \return Numerical partial of the torque w.r.t. body rotational state elements
 */
Eigen::MatrixXd calculateTorqueWrtRotationalStatePartials(
        std::function< void( Eigen::Vector7d ) > setBodyRotationalState,
        std::shared_ptr< basic_astrodynamics::TorqueModel > torqueModel,
        Eigen::Vector7d originalRotationalState,
        Eigen::VectorXd statePerturbations,
        int startIndex,
        int numberOfEntries,
        std::function< void( ) > updateFunction = emptyFunction,
        const double evaluationTime = TUDAT_NAN );

//! Function to numerical compute the partial derivative of a acceleration w.r.t. a body rotational quaternion.
/*!
 * Function to numerical compute the partial derivative of an acceleration w.r.t. a body rotational quaternion elements
 * using a first-order central difference method.
 * \param setBodyRotationalState Function to reset the current rotational state w.r.t. which the partial is to be computed.
 * \param accelerationModel Acceleration model for which the partial derivative is to be computed.
 * \param originalRotationalState Nominal rotational state at which the partial derivative is to be computed.
 * \param commandedQuaternionPerturbation Perturbation to the quaternion that is to be used
 * \param appliedQuaternionPerturbation Quaternion perturbations that were actually applied (after manifold correction)
 * \param updateFunction Function to update the required environment models following the change of the body state.
 * \param evaluationTime Time at which partial is to be evaluated (default NaN).
 * \return Numerical partial of the acceleration w.r.t. body rotational quaternion elements
 */
Eigen::MatrixXd calculateAccelerationDeviationDueToOrientationChange(
        const std::function< void( Eigen::Vector7d ) > setBodyRotationalState,
        std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        const Eigen::Vector7d& originalRotationalState,
        const Eigen::Vector4d& commandedQuaternionPerturbation,
        std::vector< Eigen::Vector4d >& appliedQuaternionPerturbation,
        std::function< void( ) > updateFunction = emptyFunction,
        const double evaluationTime = TUDAT_NAN );

//! Function to numerical compute the partial derivative of a torque w.r.t. a body translational state.
/*!
 * Function to numerical compute the partial derivative of an torque w.r.t. a body translational state elements
 * using a first-order central difference method.
 * \param setBodyState Function to reset the current translational state w.r.t. which the partial is to be computed.
 * \param torqueModel Torque model for which the partial derivative is to be computed.
 * \param originalState Nominal translational state at which the partial derivative is to be computed.
 * \param statePerturbation Perturbation to the states that is to be used
 * \param startIndex Start index in the state vector from where the statePerturbation is to be added (i.e. 0 if
 * function is to compute partial w.r.t. position, 3 if w.r.t velocity).
 * \param updateFunction Function to update the required environment models following the change of the body state.
 * \param evaluationTime Time at which partial is to be evaluated (default NaN).
 * \return Numerical partial of the torque w.r.t. body translational state elements
 */
Eigen::MatrixXd calculateTorqueWrtTranslationalStatePartials(
        std::function< void( Eigen::Vector6d ) > setBodyState,
        std::shared_ptr< basic_astrodynamics::TorqueModel > torqueModel,
        Eigen::Vector6d originalState,
        Eigen::Vector3d statePerturbation,
        int startIndex,
        std::function< void( ) > updateFunction = emptyFunction,
        const double evaluationTime = TUDAT_NAN );

//! Function to numerical compute the partial derivative of a torque w.r.t. a body rotational quaternion.
/*!
 * Function to numerical compute the partial derivative of an torque w.r.t. a body rotational quaternion elements
 * using a first-order central difference method.
 * \param setBodyRotationalState Function to reset the current rotational state w.r.t. which the partial is to be computed.
 * \param torqueModel Torque model for which the partial derivative is to be computed.
 * \param originalRotationalState Nominal rotational state at which the partial derivative is to be computed.
 * \param commandedQuaternionPerturbation Perturbation to the quaternion that is to be used
 * \param appliedQuaternionPerturbation Quaternion perturbations that were actually applied (after manifold correction)
 * \param updateFunction Function to update the required environment models following the change of the body state.
 * \param evaluationTime Time at which partial is to be evaluated (default NaN).
 * \return Numerical partial of the torque w.r.t. body rotational quaternion elements
 */
Eigen::MatrixXd calculateTorqueDeviationDueToOrientationChange(
        const std::function< void( Eigen::Vector7d ) > setBodyRotationalState,
        const std::shared_ptr< basic_astrodynamics::TorqueModel > torqueModel,
        const Eigen::Vector7d& originalRotationalState,
        const Eigen::Vector4d& commandedQuaternionPerturbation,
        std::vector< Eigen::Vector4d >& appliedQuaternionPerturbation,
        std::function< void( ) > updateFunction = emptyFunction,
        const double evaluationTime = TUDAT_NAN );

//! Function to numerical compute the partial derivative of an acceleration w.r.t. a double parameter
/*!
 * Function to numerical compute the partial derivative of an acceleration w.r.t. a double parameter,
 * using a first-order central difference method.
 * \param parameter Object describing the parameter w.r.t. which the partial is to be taken.
 * \param accelerationModel Acceleration model for which the partial derivative is to be computed.
 * \param parameterPerturbation Perturbation to be used for parameter value.
 * \param updateDependentVariables  Function to update the required environment models following the change in parameter,
 * for models that do not explicitly depend on the current time.
 * \param currentTime Time at which partial is to be computed.
 * \param timeDependentUpdateDependentVariables Function to update the required environment models following the change in
 * parameters for models that do  explicitly depend on the current time.
 * \return Numerical partial of the acceleration w.r.t. given parameter.
 */
Eigen::Vector3d calculateAccelerationWrtParameterPartials(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
        std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        double parameterPerturbation,
        std::function< void( ) > updateDependentVariables = emptyFunction,
        const double currentTime = 0.0,
        std::function< void( const double ) > timeDependentUpdateDependentVariables = emptyTimeFunction );

//! Function to numerical compute the partial derivative of an torque w.r.t. a double parameter
/*!
 * Function to numerical compute the partial derivative of an torque w.r.t. a double parameter,
 * using a first-order central difference method.
 * \param parameter Object describing the parameter w.r.t. which the partial is to be taken.
 * \param torqueModel Torque model for which the partial derivative is to be computed.
 * \param parameterPerturbation Perturbation to be used for parameter value.
 * \param updateDependentVariables  Function to update the required environment models following the change in parameter,
 * for models that do not explicitly depend on the current time.
 * \param currentTime Time at which partial is to be computed.
 * \param timeDependentUpdateDependentVariables Function to update the required environment models following the change in
 * parameters for models that do  explicitly depend on the current time.
 * \return Numerical partial of the torque w.r.t. given parameter.
 */
Eigen::Vector3d calculateTorqueWrtParameterPartials(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
        std::shared_ptr< basic_astrodynamics::TorqueModel > torqueModel,
        double parameterPerturbation,
        std::function< void( ) > updateDependentVariables = emptyFunction,
        const double currentTime = 0.0,
        std::function< void( const double ) > timeDependentUpdateDependentVariables = emptyTimeFunction );

double calculateMassRateWrtParameterPartials(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
        std::shared_ptr< basic_astrodynamics::MassRateModel > massRateModel,
        double parameterPerturbation,
        std::function< void( ) > updateDependentVariables = emptyFunction,
        const double currentTime = 0.0,
        std::function< void( const double ) > timeDependentUpdateDependentVariables = emptyTimeFunction );

//! Function to numerical compute the partial derivative of an acceleration w.r.t. a vector parameter
/*!
 * Function to numerical compute the partial derivative of an acceleration w.r.t. a vector parameter,
 * using a first-order central difference method.
 * \param parameter Object describing the parameter w.r.t. which the partial is to be taken.
 * \param accelerationModel Acceleration model for which the partial derivative is to be computed.
 * \param parameterPerturbation Perturbations to be used for parameter value.
 * \param updateDependentVariables  Function to update the required environment models following the change in parameter,
 * for models that do not explicitly depend on the current time.
 * \param currentTime Time at which partial is to be computed.
 * \param timeDependentUpdateDependentVariables Function to update the required environment models following the change in
 * parameters for models that do  explicitly depend on the current time.
 * \return Numerical partial of the acceleration w.r.t. given parameter.
 */
Eigen::Matrix< double, 3, Eigen::Dynamic > calculateAccelerationWrtParameterPartials(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
        std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        Eigen::VectorXd parameterPerturbation,
        std::function< void( ) > updateDependentVariables = emptyFunction,
        const double currentTime = 0.0,
        std::function< void( const double ) > timeDependentUpdateDependentVariables = emptyTimeFunction );

//! Function to numerical compute the partial derivative of an torque w.r.t. a vector parameter
/*!
 * Function to numerical compute the partial derivative of an torque w.r.t. a vector parameter,
 * using a first-order central difference method.
 * \param parameter Object describing the parameter w.r.t. which the partial is to be taken.
 * \param torqueModel Torque model for which the partial derivative is to be computed.
 * \param parameterPerturbation Perturbation to be used for parameter value.
 * \param updateDependentVariables  Function to update the required environment models following the change in parameter,
 * for models that do not explicitly depend on the current time.
 * \param currentTime Time at which partial is to be computed.
 * \param timeDependentUpdateDependentVariables Function to update the required environment models following the change in
 * parameters for models that do  explicitly depend on the current time.
 * \return Numerical partial of the torque w.r.t. given parameter.
 */
Eigen::Matrix< double, 3, Eigen::Dynamic > calculateTorqueWrtParameterPartials(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
        std::shared_ptr< basic_astrodynamics::TorqueModel > torqueModel,
        Eigen::VectorXd parameterPerturbation,
        std::function< void( ) > updateDependentVariables = emptyFunction,
        const double currentTime = 0.0,
        std::function< void( const double ) > timeDependentUpdateDependentVariables = emptyTimeFunction );

} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_NUMERICALACCELERATIONPARTIAL_H

/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include <boost/function.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Basics/basicTypedefs.h"

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
        boost::function< void( Eigen::Vector6d ) > setBodyState,
        boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        Eigen::Vector6d originalState,
        Eigen::Vector3d statePerturbation,
        int startIndex,
        boost::function< void( ) > updateFunction = emptyFunction,
        const double evaluationTime = TUDAT_NAN );

//! Function to numerical compute the partial derivative of an acceleration w.r.t. a double parameter
/*!
 * Function to numerical compute the partial derivative of an acceleration w.r.t. a double parameter,
 * using a first-order central difference method.
 * \param parameter Object describing the parameter w.r.t. which teh partial is to be taken.
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
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
        boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        double parameterPerturbation,
        boost::function< void( ) > updateDependentVariables = emptyFunction,
        const double currentTime = 0.0,
        boost::function< void( const double ) > timeDependentUpdateDependentVariables = emptyTimeFunction );

//! Function to numerical compute the partial derivative of an acceleration w.r.t. a vector parameter
/*!
 * Function to numerical compute the partial derivative of an acceleration w.r.t. a vector parameter,
 * using a first-order central difference method.
 * \param parameter Object describing the parameter w.r.t. which teh partial is to be taken.
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
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
        boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        Eigen::VectorXd parameterPerturbation,
        boost::function< void( ) > updateDependentVariables = emptyFunction,
        const double currentTime = 0.0,
        boost::function< void( const double ) > timeDependentUpdateDependentVariables = emptyTimeFunction );


} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_NUMERICALACCELERATIONPARTIAL_H

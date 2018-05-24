/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_QUATERNION_HISTORY_MANIPULATION_H
#define TUDAT_QUATERNION_HISTORY_MANIPULATION_H

#include <cmath>
#include <limits>
#include <map>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Basics/basicTypedefs.h"

#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"

namespace tudat
{

namespace orbital_element_conversions
{

//! Transform quaternion to opposite rotation, based on discontinuities in derivatives.
/*!
 *  Transform quaternion to opposite rotation, based on discontinuities in derivatives.
 *  \param quaternionHistoryMap Map of quaternions over time, where quaternions appear as discontinuous.
 *  \return Map of quaternions over time, where continuity is restored (returned by reference).
 */
void matchQuaternionHistory( std::map< double, Eigen::Vector4d >& quaternionHistoryMap );

//! Transform quaternion in translational or rotational state to opposite rotation, based on discontinuities
//! in derivatives.
/*!
 *  Transform quaternion in translational or rotational state to opposite rotation, based on discontinuities
 *  in derivatives.
 *  \param stateHistoryMap Map of state over time, where quaternions appear as discontinuous.
 *  \param stateType Type of state to be processed (translational or rotational).
 *  \return Map of state over time, where continuity is restored (returned by reference).
 */
void matchQuaternionHistory( std::map< double, Eigen::VectorXd >& stateHistoryMap,
                             const propagators::IntegratedStateType stateType );

//! Transform quaternion to opposite rotation, based on shadow flag of another attitude representation.
/*!
 *  Transform quaternion to opposite rotation, based on shadow flag of another attitude representation.
 *  \param quaternionHistoryMap Map of quaternions over time, where quaternions appear as discontinuous.
 *  \param attitudeHistoryMap Map of other attitude representation over time, which was used to convert to quaternions.
 *  \return Map of quaternions over time, where continuity is restored (returned by reference).
 */
void matchQuaternionHistory( std::map< double, Eigen::Vector4d >& quaternionHistoryMap,
                             const std::map< double, Eigen::Vector4d >& attitudeHistoryMap );

//! Transform quaternion in translational or rotational state to opposite rotation, based on shadow
//! flag of another attitude representation.
/*!
 *  Transform quaternion in translational or rotational state to opposite rotation, based on shadow
 *  flag of another attitude representation.
 *  \param conventionalStateHistoryMap Map of conventional state over time, where quaternions appear as discontinuous.
 *  \param propagatedStateHistoryMap Map of propagated state over time, which was used to convert to quaternions.
 *  \param stateType Type of state to be processed (translational or rotational).
 *  \return Map of state over time, where continuity is restored (returned by reference).
 */
void matchQuaternionHistory( std::map< double, Eigen::VectorXd >& conventionalStateHistoryMap,
                             const std::map< double, Eigen::VectorXd >& propagatedStateHistoryMap,
                             const propagators::IntegratedStateType stateType );

} // namespace orbital_element_conversions

} // namespace tudat

#endif // TUDAT_QUATERNION_HISTORY_MANIPULATION_H

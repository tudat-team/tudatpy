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

//! Transform quaternion to opposite rotation.
/*!
 *  Transform quaternion to opposite rotation.
 *  \param quaternionHistory Map of quaternions over time, where quaternions appear as discontinuous.
 *  \return Map of quaternions over time, where continuity is restored.
 */
void convertQuaternionHistoryToMatchSigns( std::map< double, Eigen::Vector4d >& quaternionHistoryMap );

//! Transform quaternion in translational or rotational state to opposite rotation.
/*!
 *  Transform quaternion in translational or rotational state to opposite rotation.
 *  \param quaternionHistory Map of state over time, where quaternions appear as discontinuous.
 *  \return Map of state over time, where continuity is restored.
 */
void convertQuaternionHistoryToMatchSigns( std::map< double, Eigen::VectorXd >& stateHistoryMap,
                                           const propagators::IntegratedStateType stateType );

} // namespace orbital_element_conversions

} // namespace tudat

#endif // TUDAT_QUATERNION_HISTORY_MANIPULATION_H

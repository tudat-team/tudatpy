/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ATTITUDE_ELEMENT_CONVERSIONS_H
#define TUDAT_ATTITUDE_ELEMENT_CONVERSIONS_H

#include <cmath>
#include <limits>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{

namespace orbital_element_conversions
{

//! Convert quaternions to modified Rodrigues parameters.
/*!
 *  Convert quaternions to modified Rodrigues parameters. The conversion is slightly different for modified Rodrigues parameters (MRP)
 *  than for shadow modified Rodrigues parameters (SMRP). This difference is embodied by conversionSign.
 *  \param quaternionElements Vector of quaternion elements.
 *  \return convertedModifiedRodriguesParameterElements Vector of (shadow) modified Rodrigues parameter elements.
 */
Eigen::Vector4d convertQuaternionsToModifiedRodriguesParameterElements( const Eigen::Vector4d& quaternionElements );

//! Convert modified Rodrigues parameters to quaternions.
/*!
 *  Convert modified Rodrigues parameters to quaternions. The conversion is slightly different for modified Rodrigues parameters (MRP)
 *  than for shadow modified Rodrigues parameters (SMRP). This difference is embodied by conversionSign.
 *  \param modifiedRodriguesParameterElements Vector of (shadow) modified Rodrigues parameters elements.
 *  \return convertedQuaternionElements Vector of quaternion elements.
 */
Eigen::Vector4d convertModifiedRodriguesParametersToQuaternionElements( const Eigen::Vector4d& modifiedRodriguesParameterElements );

//! Convert quaternions to exponential map.
/*!
 *  Convert quaternions to exponential map. The conversion is the same for both exponential map (EM) and shadow
 *  exponential map (SEM).
 *  \param quaternionElements Vector of quaternion elements.
 *  \return convertedExponentialMapElements Vector of (shadow) exponential map elements.
 */
Eigen::Vector4d convertQuaternionsToExponentialMapElements( const Eigen::Vector4d& quaternionElements );

//! Convert exponential map to quaternions.
/*!
 *  Convert exponential map to quaternions. The conversion is the same for both exponential map (EM) and shadow
 *  exponential map (SEM).
 *  \param exponentialMapElements Vector of (shadow) exponential map elements.
 *  \return convertedQuaternionElements Vector of quaternion elements.
 */
Eigen::Vector4d convertExponentialMapToQuaternionElements( const Eigen::Vector4d& exponentialMapElements );

} // namespace orbital_element_conversions

} // namespace tudat

#endif // TUDAT_ATTITUDE_ELEMENT_CONVERSIONS_H

/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_STELLAR_ABERRATION_CORRECTION_H
#define TUDAT_STELLAR_ABERRATION_CORRECTION_H

#include <Eigen/Core>

namespace tudat
{
namespace observation_models
{

//! Convert an astrometric line of sight to its aberrated apparent direction.
/*!
 * Applies the special-relativistic stellar-aberration transformation for an observer with the supplied inertial velocity.
 * Both vectors must be expressed in the same inertial frame.
 * \param astrometricDirection Geometric, light-time-corrected direction from observer to source. The vector need not be normalized.
 * \param observerVelocity Inertial observer velocity in m/s.
 * \return Unit apparent direction from observer to source.
 */
Eigen::Vector3d calculateApparentDirectionWithStellarAberration( const Eigen::Vector3d& astrometricDirection,
                                                                 const Eigen::Vector3d& observerVelocity );

}  // namespace observation_models
}  // namespace tudat

#endif  // TUDAT_STELLAR_ABERRATION_CORRECTION_H

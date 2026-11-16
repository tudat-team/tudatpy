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
 * Applies classical (first-order in v/c) stellar aberration for an observer with the supplied inertial velocity: the returned unit
 * vector a is the unit vector for which a - v/c is parallel to the astrometric direction u, i.e. the inverse of the aberration
 * correction u ~ a - v/c. Writing beta = v/c, this gives a = s u + beta with s = -u.beta + sqrt( 1 - beta.beta + ( u.beta )^2 ).
 * This differs from the special-relativistic transformation (e.g. Klioner, 2003, AJ 125, 1580, doi:10.1086/367593, as implemented
 * in SOFA iauAb) at order (v/c)^2, i.e. by up to about 2.5E-9 rad (0.5 mas) for an observer moving with Earth's orbital velocity.
 * Both vectors must be expressed in the same inertial frame.
 * \param astrometricDirection Geometric, light-time-corrected direction from observer to source. The vector need not be normalized.
 * \param observerVelocity Inertial observer velocity in m/s.
 * \return Unit apparent direction from observer to source.
 */
Eigen::Vector3d calculateApparentDirectionWithStellarAberration( const Eigen::Vector3d& astrometricDirection,
                                                                 const Eigen::Vector3d& observerVelocity );

//! Calculate the partial derivatives of the aberrated apparent direction.
/*!
 * Computes the Jacobians of the unit vector returned by calculateApparentDirectionWithStellarAberration with respect to its
 * (unnormalized) astrometric direction and the observer velocity, for use in observation partials.
 * \param astrometricDirection Geometric, light-time-corrected direction from observer to source. The vector need not be normalized.
 * \param observerVelocity Inertial observer velocity in m/s.
 * \param partialWrtAstrometricDirection Partial of the apparent unit direction w.r.t. astrometricDirection (returned by reference).
 * \param partialWrtObserverVelocity Partial of the apparent unit direction w.r.t. observerVelocity, in s/m (returned by reference).
 */
void calculateStellarAberrationPartials( const Eigen::Vector3d& astrometricDirection,
                                         const Eigen::Vector3d& observerVelocity,
                                         Eigen::Matrix3d& partialWrtAstrometricDirection,
                                         Eigen::Matrix3d& partialWrtObserverVelocity );

}  // namespace observation_models
}  // namespace tudat

#endif  // TUDAT_STELLAR_ABERRATION_CORRECTION_H

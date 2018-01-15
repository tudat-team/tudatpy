/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_LORENTZ_STATIC_MAGNETIC_FORCE_H
#define TUDAT_LORENTZ_STATIC_MAGNETIC_FORCE_H

#include <Eigen/Core>

namespace tudat
{
namespace electro_magnetism
{

//! Compute Lorentz force due to static magnetic field
/*!
 * Computes Lorentz force due to static magnetic field on a particle in an
 * inertial reference frame.
 * Assumes the particle is point-like.The following equation is used to calculate the force:
 * \f[
 *                \bar{f} = q * \bar{v} x \bar{B}
 * \f] 
 * where \f$f\f$ is the force on the particle, \f$q\f$ is the charge of the accelerating
 * particle, \f$\bar{v}\f$ is the velocity of the accelerating particle, and \f$\bar{B}\f$
 * is the local magnetic field. Note: This force does not take into account
 * the Coulomb force between the particle and the source.
 * \param velocityOfBodySubjectToAcceleration Velocity of body which is being accelerated by the Lorentz force. [m/s]
 * \param chargeOfBodySubjectToAcceleration Charge of body which is being accelerated by Lorentz force.         [C]
 * \param localMagneticField local magnetic field at position of body subject to acceleration                   [Tm]
 * \return Lorentz force due to static magnetic field                                                           [N]
 */

Eigen::Vector3d computeLorentzForceDueToStaticMagneticField(
        const Eigen::Vector3d& velocityOfBodySubjectToAcceleration,
        const Eigen::Vector3d& localMagneticField,
        const double chargeOfBodySubjectToAcceleration );

} // namespace electro_magnetism
} // namespace tudat

#endif // TUDAT_LORENTZ_STATIC_MAGNETIC_FORCE_H

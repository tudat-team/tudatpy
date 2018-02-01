/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *   
 *
 */

#ifndef TUDAT_LORENTZ_STATIC_MAGNETIC_ACCELERATION_H
#define TUDAT_LORENTZ_STATIC_MAGNETIC_ACCELERATION_H

#include <Eigen/Core>

namespace tudat
{
namespace electro_magnetism
{

//! Compute Lorentz acceleration due to static magnetic field
/*!
 * Computes Lorentz acceleration due to static magnetic field on a particle in an
 * inertial reference frame.
 * Assumes the particle is point-like. Note: This function does not take into account any
 * the Coulomb force between the particle and the source.
 * \param velocityOfBodySubjectToAcceleration Velocity of body which is being accelerated by the Lorentz force.  [m/s]
 * \param localMagneticField local magnetic field at position of body subject to acceleration                    [Tm]
 * \param chargeOfBodySubjectToAcceleration Charge of body which is being accelerated by Lorentz force.          [C]
 * \param massOfBodySubjectToAcceleration Mass of body which is being accelerated by Lorentz force               [kg]
 * \return Lorentz acceleration due to static magnetic field.                                                    [N]
 * \sa computeLorentzForceDueToStaticMagneticField.
 */
Eigen::Vector3d computeLorentzAccelerationDueToStaticMagneticField(
        const Eigen::Vector3d& velocityOfBodySubjectToAcceleration,
        const Eigen::Vector3d& localMagneticField,
        const double chargeOfBodySubjectToAcceleration,
        const double massOfBodySubjectToAcceleration );

} // namespace electro_magnetism
} // namespace tudat

#endif // TUDAT_LORENTZ_STATIC_MAGNETIC_ACCELERATION_H

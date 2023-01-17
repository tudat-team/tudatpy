/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */
//RP-OLD

#ifndef TUDAT_CANNON_BALL_RADIATION_PRESSURE_FORCE_H
#define TUDAT_CANNON_BALL_RADIATION_PRESSURE_FORCE_H

#include <Eigen/Core>

namespace tudat
{
namespace electromagnetism
{

//! Compute radiation pressure force using a cannon-ball model.
/*!
 * Computes radiation pressure force using a cannon-ball model, i.e., assuming force to be in
 * opposite direction of the vector to the source.
 * \param radiationPressure Radiation pressure at target. N.B: the usual way of computing the
 *          radiation pressure at the target, in case the source is the Sun, is to take the
 *          radiation presssure at 1 AU and scale it based on the distance from the Sun     [N/m^2]
 * \param vectorToSource Vector pointing from target to source. N.B: this must be a unit
 *          vector! To compute the unit vector based on a given position vector, you can
 *          use the .normalize() or .normalized() member functions of an Eigen::Vector3d
 *          object.                                                                             [-]
 * \param area Area on which radiation pressure is assumed to act.                            [m^2]
 * \param radiationPressureCoefficient Coefficient to scale effective force. Equal to 1 +
 *          emissivitty, assuming no diffuse reflection.                                        [-]
 * \return Force due to radiation pressure.                                                     [N]
 */
Eigen::Vector3d computeCannonBallRadiationPressureForce(
        const double radiationPressure,
        const Eigen::Vector3d& vectorToSource,
        const double area,
        const double radiationPressureCoefficient );

} // namespace electromagnetism
} // namespace tudat

#endif // TUDAT_CANNON_BALL_RADIATION_PRESSURE_FORCE_H

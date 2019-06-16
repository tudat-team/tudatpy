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

#ifndef TUDAT_SOLAR_SAIL_FORCE_H
#define TUDAT_SOLAR_SAIL_FORCE_H

#include <Eigen/Core>

namespace tudat
{
namespace electro_magnetism
{

//! Compute solar sail force using a non-ideal reflective model.
/*!
 * Computes solar sail force using a non-ideal reflective model
 * \param frontEmissivity Parameter determining the emissivity of the front of the sail         [-]
 * \param backEmissivity  Parameter determining the emissivity of the back of the sail          [-]
 * \param frontLambertianCoefficient  Parameter determining the Lambertian coefficient
 *        of the front of the sail                                                              [-]
 * \param backLambertianCoefficient   Parameter determining the Lambertian coefficient
 *        of the back of the sail
 * \param reflectivityCoefficient   Coefficient determining the front reflectivity of
 *        the sail
 * \param specularReflection Coefficient  Coefficient of specular reflection                    [-]
 * \param vectorToSource Vector pointing from target to source. N.B: this must be a unit
 *          vector! To compute the unit vector based on a given position vector, you can
 *          use the .normalize() or .normalized() member functions of an Eigen::Vector3d
 *          object.
 * \param radiationPressure Radiation pressure at targert                                   [N/m^2]
 * \param area Area on which radiation pressure is assumed to act.                            [m^2]
 * \param coneAngle Sail cone angle                                                           [rad]
 * \param mass Mass of accelerated body.                                                       [kg]
 * \return Solar sail force.                                                                    [N]
 */
Eigen::Vector3d computeSolarSailForce(
        const double frontEmissivityCoefficient,
        const double backEmissivityCoefficient,
        const double frontLambertianCoefficient,
        const double backLambertianCoefficient,
        const double reflectivityCoefficient,
        const double specularReflectionCoefficient,
        const Eigen::Vector3d& normalisedVectorToSource,
        const Eigen::Vector3d& normalisedVelocityVector,
        const double radiationPressure,
        const double area,
        const double coneAngle,
        const double clockAngle);


} // namespace electro_magnetism
} // namespace tudat

#endif // TUDAT_SOLAR_SAIL_FORCE_H
/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SPHERICALSTATECONVERSIONS_H
#define TUDAT_SPHERICALSTATECONVERSIONS_H

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

namespace tudat
{

namespace orbital_element_conversions
{

//! Calculate current heading angle.
/*!
 * Calculate heading angle from velocity in vertical (LVLH) frame.
 * NOTE: This function can be used for both ground- and airspeed-based variables, with both the input and output
 * always w.r.t. the same velocity vector.
 * \param velocityInVerticalFrame Current Cartesian velocity in vertical frame.
 * \return Current heading angle.
 */
double calculateHeadingAngle(
        const Eigen::Vector3d& velocityInVerticalFrame );

//! Calculate current flight path angle. Angle is defined positive upwards.
/*!
 *  Calculate flight path angle from velocity in vertical (LVLH) frame.
 *  Angle is defined positive upwards.
 * NOTE: This function can be used for both ground- and airspeed-based variables, with both the input and output
 * always w.r.t. the same velocity vector.
 *  \param velocityInVerticalFrame Current Cartesian velocity in vertical frame.
 *  \return Current flight path angle.
 */
double calculateFlightPathAngle(
        const Eigen::Vector3d& velocityInVerticalFrame );

//! Function to convert a Cartesian to a spherical orbital state.
/*!
 * Function to convert a Cartesian to a spherical orbital state. We define the spherical orbital position by the
 * radius, latitude and longitude. The spherical orbital velocity is defined by the speed, flight path angle (positive
 * upwards) and the heading angle. Note that this is distinct from a 'mathematical' spherical state, where the velocity
 * is denoted by the radius rate, latitude rate and longitude rate.
 * The spherical orbital state is often used to define entry (initial) conditions. The order of the entries in the
 * return vector are defined by the SphericalOrbitalStateElementIndices enum.
 * NOTE: This function can be used for both ground- and airspeed-based variables, with both the input and output
 * always w.r.t. the same velocity vector.
 * \param bodyFixedCartesianState Vehicle state in a frame fixed to the body w.r.t. which the state is to be computed.
 * \return Spherical orbital state representation of bodyFixedCartesianState
 */
Eigen::Vector6d convertCartesianToSphericalOrbitalState(
        const Eigen::Vector6d& bodyFixedCartesianState );

//! Function to convert a spherical orbital to a Cartesian state.
/*!
 * Function to convert a spherical orbital to a Cartesian state. We define the spherical orbital position by the
 * radius, latitude and longitude. The spherical orbital velocity is defined by the speed, flight path angle (positive
 * upwards) and the heading angle. Note that this is distinct from a 'mathematical' spherical state, where the velocity
 * is denoted by the radius rate, latitude rate and longitude rate.
 * The spherical orbital state is often used to define entry (initial) conditions. The order of the entries in the
 * return vector are defined by the SphericalOrbitalStateElementIndices enum.
 * \param sphericalOrbitalState Vehicle spherical orbital state in a frame fixed to the body w.r.t. which the state is to be
 * computed.
 * NOTE: This function can be used for both ground- and airspeed-based variables, with both the input and output
 * always w.r.t. the same velocity vector.
 * \return Cartesian state representation of sphericalOrbitalState (in same frame).
 */
Eigen::Vector6d convertSphericalOrbitalToCartesianState(
        const Eigen::Vector6d& sphericalOrbitalState );

} // namespace orbital_element_conversions

} // namespace tudat


#endif // TUDAT_SPHERICALSTATECONVERSIONS_H

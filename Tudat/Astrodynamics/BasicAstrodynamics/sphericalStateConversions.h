/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SPHERICALSTATECONVERSIONS_CPP
#define TUDAT_SPHERICALSTATECONVERSIONS_CPP

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"
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
 * \param velocityInVerticalFrame Current Cartesian velocity in vertical frame.
 * \return Current heading angle.
 */
double calculateHeadingAngle( const Eigen::Vector3d& velocityInVerticalFrame );

//! Calculate current flight path angle. Angle is defined positive upwards.
/*!
 *  Calculate flight path angle from velocity in vertical (LVLH) frame.
 *  Angle is defined positive upwards.
 *  \param velocityInVerticalFrame Current Cartesian velocity in vertical frame.
 *  \return Current flight path angle.
 */
double calculateFlightPathAngle( const Eigen::Vector3d& velocityInVerticalFrame );

basic_mathematics::Vector6d convertCartesianToSphericalOrbitalState(
        const basic_mathematics::Vector6d& bodyFixedCartesianState );

basic_mathematics::Vector6d convertSphericalOrbitalToCartesianState(
        const basic_mathematics::Vector6d& sphericalOrbitalState );
}

}
#endif // TUDAT_SPHERICALSTATECONVERSIONS_CPP

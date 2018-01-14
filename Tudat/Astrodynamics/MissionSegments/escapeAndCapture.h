/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      Note that it is only possible to define the parking orbit using the semi-major axis and
 *      eccentricity right now. One may want to add more options in the future.
 *
 */

#ifndef TUDAT_ESCAPE_AND_CAPTURE_H
#define TUDAT_ESCAPE_AND_CAPTURE_H

namespace tudat
{
namespace mission_segments
{

//! Compute escape or capture deltaV budget.
/*!
 * Calculates the required deltaV to perform a certain escape or capture maneuver to/from the
 * specified orbit from/to the specified excessVelocity.
 * \param gravitationalParameter Gravitational parameter of the escape/capture body.     [m^3 s^-2]
 * \param semiMajorAxis Semi major axis of the parking orbit.                                   [m]
 * \param eccentricity Eccentricity of the parking orbit.                                       [-]
 * \param excessVelocity Excess velocity after the escape or before the capture maneuver.  [m s^-1]
 * \return deltaV The deltaV required for the escape or capture maneuver.                  [m s^-1]
 */
double computeEscapeOrCaptureDeltaV( const double gravitationalParameter,
                                     const double semiMajorAxis,
                                     const double eccentricity,
                                     const double excessVelocity );

} // namespace mission_segments
} // namespace tudat

#endif // TUDAT_ESCAPE_AND_CAPTURE_H

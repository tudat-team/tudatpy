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
 *      Wakker, K. F. (2007), Lecture Notes Astrodynamics II (Chapter 18), TU Delft course AE4-874,
 *          Delft University of technology, Delft, The Netherlands.
 *
 */

#include <cmath>

#include "Tudat/Astrodynamics/MissionSegments/escapeAndCapture.h"

namespace tudat
{
namespace mission_segments
{
//! Compute escape or capture deltaV budget.
double computeEscapeOrCaptureDeltaV( const double gravitationalParameter,
                                     const double semiMajorAxis,
                                     const double eccentricity,
                                     const double excessVelocity )
{
    // Calculate the pericenter radius from the semi-major axis and eccentricity.
    const double pericenterRadius = semiMajorAxis * ( 1.0 - eccentricity );

    // Calculate deltaV using Equation 18-28 of [Wakker, 2007].
    return std::sqrt( 2.0 * gravitationalParameter / pericenterRadius +
                      excessVelocity * excessVelocity ) -
           std::sqrt( 2.0 * gravitationalParameter / pericenterRadius -
                      gravitationalParameter / semiMajorAxis );
}

} // namespace mission_segments
} // namespace tudat

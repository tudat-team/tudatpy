/*    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110129    E. Iorfida        Creation of code.
 *      110131    E. Iorfida        Added comments and pointerToCelestialBody.
 *      110202    E. Iorfida        Modified structure of the code, unique base class for launch
 *                                  and capture paths.
 *      110206    E. Iorfida        Modified some comments and name of base class to
 *                                  EscapeAndCapture.
 *      110208    E. Iorfida        Deleted inheritance from TrajectoryDesignMethod, and
 *                                  execute( ), function too. Modified getDeltaV into
 *                                  computeDeltaV.
 *      110214    E. Iorfida        Deleted temporary centralBodyRadius, replaced by an element of
 *                                  GeometricShapes.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120531    P. Musegaas       Code completely rewritten. Made it a free function.
 *      120625    P. Musegaas       Minor changes.
 *
 *    References
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

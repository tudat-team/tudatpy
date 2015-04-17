/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      100906    K. Kumar          First creation of code.
 *      111115    K. Kumar          Added checker info; corrected Doxygenc comments.
 *      120127    D. Dirkx          Moved file to Tudat core.
 *      121205    D. Dirkx          Migrated namespace to directory-based protocol and added
 *                                  backwards compatibility.
 *
 *    References
 *
 *    Notes
 *      Backwards compatibility of namespaces is implemented for Tudat Core 2 in this file. The
 *      code block marked "DEPRECATED!" at the end of the file should be removed in Tudat Core 3.
 *
 */

#ifndef TUDAT_ASTRODYNAMICS_FUNCTIONS_H
#define TUDAT_ASTRODYNAMICS_FUNCTIONS_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

namespace tudat
{
namespace basic_astrodynamics
{

//! Compute two-body orbital period.
/*!
 * Computes the two-body orbital period of an orbiting body that follows a closed conic section
 * (circle or ellipse Kepler orbit). The mass of the orbiting body is set to that of a test
 * particle by default.
 * \param semiMajorAxis Semi-major axis of Kepler orbit (circle or ellipse).
 * \param gravitationalParameterOfCentralBody Gravitational parameter of central body.
 * \param massOfOrbitingBody Mass of orbiting body.
 * \return Two-body orbital period.
 */
double computeKeplerOrbitalPeriod( const double semiMajorAxis,
                                   const double gravitationalParameterOfCentralBody,
                                   const double massOfOrbitingBody = 0.0 );

//! Compute two-body angular momentum.
/*!
 * Computes the angular momentum of an orbiting body that follows a conic section (Kepler orbit), 
 * relative to the center-of-mass of the central body. The default mass value is for the angular 
 * momentum per unit mass.
 * \param semiMajorAxis Semi-major axis of Kepler orbit.
 * \param eccentricity Eccentricity of Kepler orbit.
 * \param gravitationalParameterOfCentralBody Gravitational parameter of central body.
 * \param massOfOrbitingBody Mass of orbiting body.
 * \return Two-body angular momentum.
 */
double computeKeplerAngularMomentum( const double semiMajorAxis, const double eccentricity,
                                     const double gravitationalParameterOfCentralBody,
                                     const double massOfOrbitingBody = 1.0 );

//! Compute two-body mean motion.
/*!
 * Computes the two-body mean motion of an orbiting body that follows a conic section
 * (Kepler orbit). The mass of the orbiting body is set to that of a test particle by default.
 * \param semiMajorAxis Semi-major axis of Kepler orbit.
 * \param gravitationalParameterOfCentralBody Gravitational parameter of central body.
 * \param massOfOrbitingBody Mass of orbiting body.
 * \return Two-body mean motion.
 */
double computeKeplerMeanMotion( const double semiMajorAxis,
                                 const double gravitationalParameterOfCentralBody,
                                 const double massOfOrbitingBody = 0.0 );

//! Compute Kepler energy.
/*!
 * Computes the energy of an orbiting body that follows a conic section (Kepler orbit). The 
 * default mass value is for the two-body orbital energy per unit mass. For closed conic sections 
 * (circles, ellipses), the semi-major axis is positive, and for open sections (hyperbolas) the 
 * semi-major axis is negative.
 * \param semiMajorAxis Semi-major axis of Kepler orbit.
 * \param gravitationalParameterOfCentralBody Gravitational parameter of central body.
 * \param massOfOrbitingBody Mass of orbiting body.
 * \return Kepler orbital energy.
 */
double computeKeplerEnergy( const double semiMajorAxis,
                            const double gravitationalParameterOfCentralBody,
                            const double massOfOrbitingBody = 1.0 );

//! Compute synodic period.
/*!
 * Computes synodic period between two bodies in different Kepler orbits (closed conic sections).
 * The orbital periods must be positive values for the synodic period to be sensible.
 * \param orbitalPeriodBody1 Orbital period of Body 1.
 * \param orbitalPeriodBody2 Orbital period of Body 2.
 * \return Synodic period.
 */
double computeSynodicPeriod( const double orbitalPeriodBody1, const double orbitalPeriodBody2 );

} // namespace basic_astrodynamics
} // namespace tudat


#endif // TUDAT_ASTRODYNAMICS_FUNCTIONS_H

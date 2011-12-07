/*! \file astrodynamicsFunctions.h
 *    This header file contains general astrodynamics functions.
 *
 *    Path              : /Astrodynamics/
 *    Version           : 2
 *    Check status      : Unchecked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 11 November, 2011
 *    Last modified     : 15 November, 2011
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      100906    K. Kumar          First creation of code.
 *      111115    K. Kumar          Added checker info; corrected Doxygenc comments.
 */

#ifndef ASTRODYNAMICSFUNCTIONS_H
#define ASTRODYNAMICSFUNCTIONS_H

// Include statements.
#include "Astrodynamics/physicalConstants.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
{

//! Astrodynamics namespace.
/*!
 * Astrodynamics namespace.
 */
namespace astrodynamics
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
double computeTwoBodyOrbitalPeriod( double semiMajorAxis,
                                    double gravitationalParameterOfCentralBody,
                                    double massOfOrbitingBody = 0.0 );

//! Compute two-body angular momentum.
/*!
 * Computes the two-body angular momentum of an orbiting body that follows a conic section
 * (Kepler orbit). The default mass value is for the angular momentum per unit mass.
 * \param semiMajorAxis Semi-major axis of Kepler orbit.
 * \param eccentricity Eccentricity of Kepler orbit.
 * \param gravitationalParameterOfCentralBody Gravitational parameter of central body.
 * \param massOfOrbitingBody Mass of orbiting body.
 * \return Two-body angular momentum.
 */
double computeTwoBodyAngularMomentum( double semiMajorAxis, double eccentricity,
                                      double gravitationalParameterOfCentralBody,
                                      double massOfOrbitingBody = 1.0 );

//! Compute two-body mean motion.
/*!
 * Computes the two-body mean motion of an orbiting body that follows a conic section
 * (Kepler orbit). The mass of the orbiting body is set to that of a test particle by default.
 * \param semiMajorAxis Semi-major axis of Kepler orbit.
 * \param gravitationalParameterOfCentralBody Gravitational parameter of central body.
 * \param massOfOrbitingBody Mass of orbiting body.
 * \return Two-body mean motion.
 */
double computeTwoBodyMeanMotion( double semiMajorAxis, double gravitationalParameterOfCentralBody,
                                 double massOfOrbitingBody = 0.0 );

//! Compute two-body orbital energy.
/*!
 * Computes the two-body orbital energy of an orbiting body that follows a conic section
 * (Kepler orbit). The default mass value is for the two-body orbital enery per unit mass. For
 * closed conic sections (circles, ellipses), the semi-major axis is positive, and for open
 * sections (hyperbolas) the semi-major axis is negative.
 * \param semiMajorAxis Semi-major axis of Kepler orbit.
 * \param gravitationalParameterOfCentralBody Gravitational parameter of central body.
 * \param massOfOrbitingBody Mass of orbiting body.
 * \return Two-body orbital energy.
 */
double computeTwoBodyOrbitalEnergy( double semiMajorAxis,
                                    double gravitationalParameterOfCentralBody,
                                    double massOfOrbitingBody = 1.0 );

//! Compute synodic period.
/*!
 * Computes synodic period between two bodies in different Kepler orbits (closed conic sections).
 * The orbital periods must be positive values for the synodic period to be sensible.
 * \param orbitalPeriodBody1 Orbital period of Body 1.
 * \param orbitalPeriodBody2 Orbital period of Body 2.
 * \return Synodic period.
 */
double computeSynodicPeriod( double orbitalPeriodBody1, double orbitalPeriodBody2 );

}

}

#endif // ASTRODYNAMICSFUNCTIONS_H

// End of file.

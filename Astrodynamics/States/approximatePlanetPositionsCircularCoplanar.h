/*! \file approximatePlanetPositionsCircularCoplanar.h
 *    This header file contains the definition of an ephemeris class that makes use of the JPL
 *   "Approximate Positions of Major Planets" ( http://ssd.jpl.nasa.gov/?planet_pos ) to retrieve
 *    initial ephemeris data for a specific planet. The ephemeris is valid for circular, coplanar
 *    orbits.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : L. van der Ham
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : L.vanderHam@student.tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 29 June, 2011
 *    Last modified     : 3 August, 2011
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the
 *          Major Planets, http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf,
 *          last accessed: 24 February, 2011.
 *
 *    Notes
 *      The circular coplanar orbits apply the same orbital elements as the
 *      3D case, but the eccentricity and the inclination are implicitly set to zero.
 *      Also, the elements do not vary in time, except for the anomalies of course.
 *      By implicitly setting the inclination equal to zero, all planetary orbits
 *      lie in the ecliptic plane.
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
 *      110629    L. van der Ham    First creation of code.
 *      110803    L. van der Ham    Separated this code from approximatePlanetPositions.
 */

#ifndef APPROXIMATEPLANETPOSITIONSCIRCULARCOPLANAR_H
#define APPROXIMATEPLANETPOSITIONSCIRCULARCOPLANAR_H

// Include statements.
#include "Astrodynamics/States/approximatePlanetPositionsBase.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
{

//! Ephemeris class using JPL "Approximate Positions of Major Planets".
/*!
 * Ephemeris class using JPL "Approximate Positions of Major Planets".
 */

class ApproximatePlanetPositionsCircularCoplanar : public ApproximatePlanetPositionsBase
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    ApproximatePlanetPositionsCircularCoplanar( ) : constantOrbitalRadius_( -0.0 ) { }

    //! Get state from ephemeris; circular, coplanar case.
    /*!
     * Returns state in Cartesian elements from ephemeris for circular and coplanar orbit.
     * \return State in Cartesian elements from ephemeris for circular and coplanar orbit.
     */
    CartesianElements* getStateFromEphemeris( double julianDate );

protected:

private:

    //! Orbital radius.
    /*!
     * Constant orbital radius for circular orbit.
     */
    double constantOrbitalRadius_;
};

}

#endif // APPROXIMATEPLANETPOSITIONSCIRCULARCOPLANAR_H

// End of file.

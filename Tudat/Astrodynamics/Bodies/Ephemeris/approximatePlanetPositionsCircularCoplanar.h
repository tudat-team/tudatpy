/*    Copyright (c) 2010-2012 Delft University of Technology.
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
 *      120322    D. Dirkx          Modified to new Ephemeris interfaces.
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the Major Planets,
 *          http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf, last accessed: 24 February, 2011.
 *
 */

// Temporary notes (move to class/function doxygen):
// The circular coplanar orbits apply the same orbital elements as the
// 3D case, but the eccentricity and the inclination are implicitly set to zero.
// Also, the elements do not vary in time, except for the anomalies of course.
// By implicitly setting the inclination equal to zero, all planetary orbits
// lie in the ecliptic plane.
// 

#ifndef TUDAT_APPROXIMATE_PLANET_POSITIONS_CIRCULAR_COPLANAR_H
#define TUDAT_APPROXIMATE_PLANET_POSITIONS_CIRCULAR_COPLANAR_H

#include "Tudat/Astrodynamics/Bodies/Ephemeris/approximatePlanetPositionsBase.h"

namespace tudat
{

namespace ephemerides
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
    ApproximatePlanetPositionsCircularCoplanar( BodiesWithEphemerisData bodyWithEphemerisData )
        : constantOrbitalRadius_( -0.0 )
    {
        setPlanet( bodyWithEphemerisData );
    }

    //! Get state from ephemeris; circular, coplanar case.
    /*!
     * Returns state in Cartesian elements from ephemeris for circular and coplanar orbit.
     * \return State in Cartesian elements from ephemeris for circular and coplanar orbit.
     */

    Eigen::VectorXd getCartesianStateFromEphemeris( const double julianDate );

protected:

private:

    //! Orbital radius.
    /*!
     * Constant orbital radius for circular orbit.
     */
    double constantOrbitalRadius_;
};

} // namespace ephemerides
} // namespace tudat

#endif // TUDAT_APPROXIMATE_PLANET_POSITIONS_CIRCULAR_COPLANAR_H

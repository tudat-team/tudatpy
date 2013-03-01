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
 *      110629    L. van der Ham    File created.
 *      110803    L. van der Ham    Separated this code from approximatePlanetPositions.
 *      120322    D. Dirkx          Modified to new Ephemeris interfaces.
 *      130121    K. Kumar          Updated VectorXd to Vector6d; added shared-ptr typedef.
 *      130120    D. Dirkx          Updated with new Julian day + seconds since Julian day input.
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the Major Planets,
 *          http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf, last accessed: 24 February, 2011.
 *
 *    Notes
 *      The circular coplanar orbits apply the same orbital elements as the 3D case, but the
 *      eccentricity and the inclination are implicitly set to zero. Also, the elements do not vary
 *      in time, except for the anomalies of course. By implicitly setting the inclination equal to
 *      zero, all planetary orbits lie in the ecliptic plane.
 *
 */

#ifndef TUDAT_APPROXIMATE_PLANET_POSITIONS_CIRCULAR_COPLANAR_H
#define TUDAT_APPROXIMATE_PLANET_POSITIONS_CIRCULAR_COPLANAR_H

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositionsBase.h"

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

namespace tudat
{
namespace ephemerides
{

//! Ephemeris class for circular, coplanar orbits using Approximate Positions of Major Planets.
/*!
 * Ephemeris class for circular, coplanar planetary orbits, using JPL
 * "Approximate Positions of Major Planets".
 */
class ApproximatePlanetPositionsCircularCoplanar : public ApproximatePlanetPositionsBase
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    ApproximatePlanetPositionsCircularCoplanar(
            BodiesWithEphemerisData bodyWithEphemerisData,
            const double aSunGravitationalParameter = 1.32712440018e20 )
        : ApproximatePlanetPositionsBase( aSunGravitationalParameter ),
          constantOrbitalRadius_( -0.0 )
    {
        setPlanet( bodyWithEphemerisData );
    }

    //! Get state from ephemeris; circular, coplanar case.
    /*!
     * Returns state in Cartesian elements from ephemeris for circular and coplanar orbit.
     * \param secondsSinceEpoch Seconds since epoch.
     * \param julianDayAtEpoch Reference epoch in Julian day.
     * \return State in Cartesian elements from ephemeris for circular and coplanar orbit.
     */
    basic_mathematics::Vector6d getCartesianStateFromEphemeris(
            const double secondsSinceEpoch, const double julianDayAtEpoch );

protected:

private:

    //! Orbital radius.
    /*!
     * Constant orbital radius for circular orbit.
     */
    double constantOrbitalRadius_;
};

//! Typedef for shared-pointer to ApproximatePlanetPositionsCircularCoplanar object.
typedef boost::shared_ptr< ApproximatePlanetPositionsCircularCoplanar >
ApproximatePlanetPositionsCircularCoplanarPointer;

} // namespace ephemerides
} // namespace tudat

#endif // TUDAT_APPROXIMATE_PLANET_POSITIONS_CIRCULAR_COPLANAR_H

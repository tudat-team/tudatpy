/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
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

#include <memory>

#include "tudat/astro/ephemerides/approximatePlanetPositionsBase.h"

#include "tudat/basics/basicTypedefs.h"

namespace tudat
{
namespace ephemerides
{

//! Ephemeris class for circular, coplanar orbits using Approximate Positions of Major Planets.
/*!
 * Ephemeris class for circular, coplanar planetary orbits, using JPL
 * "Approximate Positions of Major Planets".
 */
class ApproximateJplCircularCoplanarEphemeris : public ApproximateJplSolarSystemEphemerisBase
{
public:

    using Ephemeris::getCartesianState;
    using Ephemeris::getCartesianLongState;


    //! Default constructor.
    /*!
     * Default constructor that initializes the class from the body for which the position is
     * approximated and the gravitational parameter of the Sun (default 1.32712440018e20). Other
     * members are initialized to default values.
     *
     * \param bodyName The body for which the position is approximated.
     * \param sunGravitationalParameter The gravitational parameter of the Sun [m^3/s^2].
     */
    ApproximateJplCircularCoplanarEphemeris(
            const std::string bodyName,
            const double sunGravitationalParameter = 1.32712440018e20 )
        : ApproximateJplSolarSystemEphemerisBase( sunGravitationalParameter ),
          constantOrbitalRadius_( -0.0 )
    {
        this->setPlanet( bodyName );
    }

    //! Get state from ephemeris; circular, coplanar case.
    /*!
     * Returns state in Cartesian elements from ephemeris for circular and coplanar orbit.
     * \param secondsSinceEpoch Seconds since epoch.
     * \return State in Cartesian elements from ephemeris for circular and coplanar orbit.
     */
    Eigen::Vector6d getCartesianState(
            const double secondsSinceEpoch );

protected:

private:

    //! Orbital radius.
    /*!
     * Constant orbital radius for circular orbit.
     */
    double constantOrbitalRadius_;

    double referenceJulianDate_;
};

//! Typedef for shared-pointer to ApproximateJplCircularCoplanarEphemeris object.
typedef std::shared_ptr< ApproximateJplCircularCoplanarEphemeris >
ApproximateSolarSystemEphemerisCircularCoplanarPointer;

} // namespace ephemerides
} // namespace tudat

#endif // TUDAT_APPROXIMATE_PLANET_POSITIONS_CIRCULAR_COPLANAR_H

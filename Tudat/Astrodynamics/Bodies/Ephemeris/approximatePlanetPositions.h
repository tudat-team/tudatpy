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
 *      110221    K. Kumar          First creation of code.
 *      110224    K. Kumar          Renamed class and file.
 *      110810    J. Leloux         Corrected doxygen documentation.
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the
 *          Major Planets, http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf,
 *          last accessed: 24 February, 2011.
 *
 */

#ifndef TUDAT_APPROXIMATE_PLANET_POSITIONS_H
#define TUDAT_APPROXIMATE_PLANET_POSITIONS_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToEccentricAnomaly.h"
#include "Tudat/Astrodynamics/Bodies/Ephemeris/approximatePlanetPositionsBase.h"
#include "Tudat/Mathematics/RootFindingMethods/newtonRaphson.h"

namespace tudat
{

//! Ephemeris class using JPL "Approximate Positions of Major Planets".
/*!
 * Ephemeris class using JPL "Approximate Positions of Major Planets".
 */
class ApproximatePlanetPositions : public ApproximatePlanetPositionsBase
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    ApproximatePlanetPositions( ) : eccentricAnomalyAtGivenJulianDate_( -0.0 ),
        longitudeOfPerihelionAtGivenJulianDate_( -0.0 ), meanAnomalyAtGivenJulianDate_( -0.0 ),
        trueAnomalyAtGivenJulianData_( -0.0 ) { }

    //! Get state from ephemeris.
    /*!
     * Returns state in Cartesian elements from ephemeris.
     * \return State in Cartesian elements from ephemeris.
     */
    CartesianElements* getStateFromEphemeris( double julianDate );

protected:

private:

    //! Eccentric anomaly at given Julian date.
    /*!
     * Eccentric anomaly of planet at given Julian date.
     */
    double eccentricAnomalyAtGivenJulianDate_;

    //! Longitude of perihelion at given Julian date.
    /*!
     * Longitude of perihelion of planet at given Julian date.
     */
    double longitudeOfPerihelionAtGivenJulianDate_;

    //! Mean anomaly at given Julian date.
    /*!
     * Mean anomaly of planet at given Julian date.
     */
    double meanAnomalyAtGivenJulianDate_;

    //! True anomaly at given Julian date.
    /*!
     * True anomaly of planet at given Julian date.
     */
    double trueAnomalyAtGivenJulianData_;

    //! Newton-Raphson method.
    /*!
     * Newton-Raphson method.
     */
    NewtonRaphson newtonRaphson_;
};

} // namespace tudat

#endif // TUDAT_APPROXIMATE_PLANET_POSITIONS_H

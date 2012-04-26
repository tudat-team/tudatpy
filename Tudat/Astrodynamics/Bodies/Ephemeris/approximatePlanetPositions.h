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
 *      120322    D. Dirkx          Modified to new Ephemeris interfaces.
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the Major Planets,
 *          http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf, last accessed: 24 February, 2011.
 *
 */

#ifndef TUDAT_APPROXIMATE_PLANET_POSITIONS_H
#define TUDAT_APPROXIMATE_PLANET_POSITIONS_H

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToEccentricAnomaly.h"
#include "Tudat/Astrodynamics/Bodies/Ephemeris/approximatePlanetPositionsBase.h"
#include "Tudat/Mathematics/RootFindingMethods/newtonRaphson.h"

namespace tudat
{
namespace ephemerides
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
    ApproximatePlanetPositions( BodiesWithEphemerisData bodyWithEphemerisData )
        : eccentricAnomalyAtGivenJulianDate_( TUDAT_NAN ),
          longitudeOfPerihelionAtGivenJulianDate_( TUDAT_NAN ),
          meanAnomalyAtGivenJulianDate_( TUDAT_NAN ),
          trueAnomalyAtGivenJulianData_( TUDAT_NAN )
    {
        setPlanet( bodyWithEphemerisData );
    }

    //! Get cartesian state from ephemeris.
    /*!
     * Returns cartesian state from ephemeris.
     * \return State in Cartesian elements from ephemeris.
     */
    Eigen::VectorXd getCartesianStateFromEphemeris( const double julianDate );

    //! Get keplerian state from ephemeris.
    /*!
     * Returns keplerian state in from ephemeris.
     * \return State in Keplerian elements from ephemeris.
     */
    Eigen::VectorXd getKeplerianStateFromEphemeris( const double julianDate );

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
};

} // namespace ephemerides
} // namespace tudat

#endif // TUDAT_APPROXIMATE_PLANET_POSITIONS_H

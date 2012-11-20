/*    Copyright (c) 2010-2012, Delft University of Technology
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
 *      110221    K. Kumar          Creation of code.
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
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositionsBase.h"

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
    ApproximatePlanetPositions( BodiesWithEphemerisData bodyWithEphemerisData,
                                const double aSunGravitationalParameter = 1.32712440018e20 )
        : ApproximatePlanetPositionsBase( aSunGravitationalParameter ),
          eccentricAnomalyAtGivenJulianDate_( TUDAT_NAN ),
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

/*! \file approximatePlanetPositions.h
 *    This header file contains the definition of an ephemeris class that makes
 *    use of the JPL "Approximate Positions of Major Planets"
 *    ( http://ssd.jpl.nasa.gov/?planet_pos ) to retrieve ephemeris data for a
 *    specific planet. The ephemeris file used is for the period 3000 BC to
 *    3000 AD.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Date created      : 21 February, 2011
 *    Last modified     : 21 February, 2011
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the
 *          Major Planets, http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf,
 *          last accessed: 24 February, 2011.
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 */

#ifndef APPROXIMATEPLANETPOSITIONS_H
#define APPROXIMATEPLANETPOSITIONS_H

// Include statements.
#include <cmath>
#include <fstream>
#include <map>
#include <string>
#include "approximatePlanetPositionsDataContainer.h"
#include "basicMathematicsFunctions.h"
#include "cartesianElements.h"
#include "convertMeanAnomalyToEccentricAnomaly.h"
#include "ephemeris.h"
#include "keplerianElements.h"
#include "newtonRaphson.h"
#include "state.h"
#include "unitConversions.h"

// Using declarations.
using std::map;

//! Ephemeris class using JPL "Approximate Positions of Major Planets".
/*!
 * Ephemeris class using JPL "Approximate Positions of Major Planets".
 */
class ApproximatePlanetPositions : public Ephemeris
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    ApproximatePlanetPositions( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~ApproximatePlanetPositions( );

    //! Set planet.
    /*!
     * Sets planet to retrieve ephemeris data for.
     * \param Planet planet.
     */
    void setPlanet( BodiesWithEphemerisData bodyWithEphemerisData );

    //! Get state from ephemeris.
    /*!
     * Returns state in Cartesian elements from ephemeris.
     * \return State in Cartesian elements from ephemeris.
     */
    CartesianElements* getStateFromEphemeris( const double& julianDate );

    //! Approximate planet positions data container.
    /*!
     * Approximate planet positions data container.
     */
    ApproximatePlanetPositionsDataContainer approximatePlanetPositionsDataContainer_;

protected:

private:

    //! Eccentric anomaly at given Julian date.
    /*!
     * Eccentric anomaly of planet at given Julian date.
     */
    double eccentricAnomalyAtGivenJulianDate_;

    //! Julian date.
    /*!
     * Julian date at which to obtain planet's orbital elements.
     */
    double julianDate_;

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

    //! Mean longitude at given Julian date.
    /*!
     * Mean longitude of planet at given Julian date.
     */
    double meanLongitudeAtGivenJulianDate_;

    //! Number of centuries passed J2000.
    /*!
     * Number of centuries passed J2000, computed using Julian date defined
     * by user when using getStateFromEphemeris( ).
     */
    double numberOfCenturiesPastJ2000_;

    //! True anomaly at given Julian date.
    /*!
     * True anomaly of planet at given Julian date.
     */
    double trueAnomalyAtGivenJulianData_;

    //! Map container of data from ephemeris file.
    /*!
     * Map container of string data from ephemeris data file.
     */
    map< unsigned int, string > containerOfDataFromEphemerisFile_;

    //! Convert mean anomaly to eccentric anomaly.
    /*!
     * Convert mean anomaly to eccentric anomaly.
     */
    orbital_element_conversions::ConvertMeanAnomalyToEccentricAnomaly
            convertMeanAnomalyToEccentricAnomaly_;

    //! Cartesian elements of planet at given Julian date.
    /*!
     * Cartesian elements of planet at given Julian date.
     */
    CartesianElements planetCartesianElementsAtGivenJulianDate_;

    //! Keplerian elements of planet at given Julian date.
    /*!
     * Keplerian elements of planet at given Julian date.
     */
    KeplerianElements planetKeplerianElementsAtGivenJulianDate_;

    //! Newton-Raphson method.
    /*!
     * Newton-Raphson method.
     */
    NewtonRaphson newtonRaphson_;

    //! Parse ephemeris line data.
    /*!
     * Parse ephemeris line data.
     * \param firstLineNumber First line number.
     */
    void parseEphemerisLineData_( const unsigned int& firstLineNumber );

    //! Parse line data for extra terms for ephemeris.
    /*!
     * Parse ephemeris line data.
     * \param lineNumber Line number.
     */
    void parseExtraTermsEphemerisLineData_( const unsigned int& lineNumber );
};

#endif // APPROXIMATEPLANETPOSITIONS_H

// End of file.

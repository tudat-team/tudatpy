/*! \file approximatePlanetPositionsBase.h
 *    This header file contains the definition of an ephemeris base class that makes use of the JPL
 *    "Approximate Positions of Major Planets" ( http://ssd.jpl.nasa.gov/?planet_pos ) to retrieve
 *    ephemeris data for a specific planet. The ephemeris file used is for the period 3000 BC to
 *    3000 AD.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 5
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : L. van der Ham
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : L.vanderHam@student.tudelft.nl
 *
 *    Checker           : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 21 February, 2011
 *    Last modified     : 24 August, 2011
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the
 *          Major Planets, http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf,
 *          last accessed: 24 February, 2011.
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
 *      110221    K. Kumar          First creation of code.
 *      110224    K. Kumar          Renamed class and file.
 *      110629    L. van der Ham    Added function coplanar circular orbits.
 *      110803    L. van der Ham    Created base class and separated approximatePlanetPositions
 *                                  from approximatePlanetPositionsCircularCoplanar.
 *      110824    J. Leloux         Corrected doxygen documentation.
 */

#ifndef APPROXIMATEPLANETPOSITIONSBASE_H
#define APPROXIMATEPLANETPOSITIONSBASE_H

// Include statements.
#include <cmath>
#include <fstream>
#include <map>
#include <string>
#include "Astrodynamics/States/approximatePlanetPositionsDataContainer.h"
#include "Astrodynamics/States/cartesianElements.h"
#include "Astrodynamics/States/ephemeris.h"
#include "Astrodynamics/States/keplerianElements.h"
#include "Mathematics/basicMathematicsFunctions.h"
#include "Mathematics/unitConversions.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
{

// Using declarations.
using std::map;

//! Ephemeris class using JPL "Approximate Positions of Major Planets".
/*!
 * Ephemeris class using JPL "Approximate Positions of Major Planets".
 */
class ApproximatePlanetPositionsBase : public Ephemeris
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    ApproximatePlanetPositionsBase( );

    //! Set planet.
    /*!
     * Sets planet to retrieve ephemeris data for.
     * \param bodyWithEphemerisData Planet.
     */
    void setPlanet( BodiesWithEphemerisData bodyWithEphemerisData );

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

    //! Get Cartesian state at the given Julian date from ephemeris.
    /*!
     * Returns state in Cartesian elements from ephemeris at the given Kulian date (UT1).
     * \param julianDate Get the ephemeris at this Julian date, UT1.
     * \return State in Cartesian elements from ephemeris.
     */
    CartesianElements* getStateFromEphemeris( const double& julianDate ) = 0;

protected:

    //! Julian date.
    /*!
     * Julian date at which to obtain planet's orbital elements.
     */
    double julianDate_;

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

    //! Map container of data from ephemeris file.
    /*!
     * Map container of string data from ephemeris data file.
     */
    map< unsigned int, string > containerOfDataFromEphemerisFile_;

    //! Approximate planet positions data container.
    /*!
     * Approximate planet positions data container.
     */
    ApproximatePlanetPositionsDataContainer approximatePlanetPositionsDataContainer_;

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

private:

};

}

#endif // APPROXIMATEPLANETPOSITIONSBASE_H

// End of file.

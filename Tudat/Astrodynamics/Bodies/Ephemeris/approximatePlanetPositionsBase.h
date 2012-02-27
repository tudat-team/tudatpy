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
 *      110629    L. van der Ham    Added function coplanar circular orbits.
 *      110803    L. van der Ham    Created base class and separated approximatePlanetPositions
 *                                  from approximatePlanetPositionsCircularCoplanar.
 *      110824    J. Leloux         Corrected doxygen documentation.
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the
 *          Major Planets, http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf,
 *          last accessed: 24 February, 2011.
 *
 */

#ifndef TUDAT_APPROXIMATE_PLANET_POSITIONS_BASE_H
#define TUDAT_APPROXIMATE_PLANET_POSITIONS_BASE_H

#include <cmath>
#include <fstream>
#include <map>
#include <string>
#include "Tudat/Astrodynamics/Bodies/Ephemeris/approximatePlanetPositionsDataContainer.h"
#include "Tudat/Astrodynamics/States/cartesianElements.h"
#include "Tudat/Astrodynamics/Bodies/Ephemeris/ephemeris.h"
#include "Tudat/Astrodynamics/States/keplerianElements.h"
#include <TudatCore/Mathematics/BasicMathematics/basicMathematicsFunctions.h>

namespace tudat
{

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
    //! Default constructor.
    ApproximatePlanetPositionsBase( )
        : julianDate_( -0.0 ),
          meanLongitudeAtGivenJulianDate_( -0.0 ),
          numberOfCenturiesPastJ2000_( -0.0 )
    { }

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
    CartesianElements* getStateFromEphemeris( double julianDate ) = 0;
    
    //! Load in ephemeris data for planets.
    /*!
     * This method opens and parses the p_elem_t2.txt ephemeris files for the planet positions.
     * The resulting data is stored in 
     * ApproximatePlanetPositionsBase::containerOfDataFromEphemerisFile_ to be used in the
     * generation of planet ephemeris. This method is automatically invoked if you call 
     * ApproximatePlanetPositionsBase::setPlanet( BodiesWithEphemerisData ).
     */
    void reloadData( );

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
    std::map< unsigned int, std::string > containerOfDataFromEphemerisFile_;

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

} // namespace tudat

#endif // TUDAT_APPROXIMATE_PLANET_POSITIONS_BASE_H

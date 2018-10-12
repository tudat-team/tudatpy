/*    Copyright (c) 2010-2018, Delft University of Technology
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
 */

#ifndef TUDAT_APPROXIMATE_PLANET_POSITIONS_BASE_H
#define TUDAT_APPROXIMATE_PLANET_POSITIONS_BASE_H

#include <cmath>
#include <map>
#include <string>

#include <memory>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositionsDataContainer.h"
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{
namespace ephemerides
{

//! Ephemeris base class using JPL "Approximate Positions of Major Planets".
/*!
 * Ephemeris base class using JPL "Approximate Positions of Major Planets".
 */
class ApproximatePlanetPositionsBase : public Ephemeris
{
public:

    //! Bodies with ephemeris data.
    enum BodiesWithEphemerisData
    {
        mercury, venus, earthMoonBarycenter, mars, jupiter, saturn, uranus, neptune, pluto
    };

    //! Default constructor.
    /*!
     * Default constructor of the base class, initializes the gravitational parameter of the Sun to
     * the input value, and all other private base class members to default values.
     *
     * \param aSunGravitationalParameter The gravitational parameter of the Sun [m^3/s^2].
     * \sa ApproximatePlanetPositions, ApproximatePlanetPositionsCircularCoplanar.
     */
    ApproximatePlanetPositionsBase( const double aSunGravitationalParameter )
        : Ephemeris( "Sun", "J2000" ),
          sunGravitationalParameter( aSunGravitationalParameter ),
          julianDate_( -0.0 ),
          meanLongitudeAtGivenJulianDate_( -0.0 ),
          numberOfCenturiesPastJ2000_( -0.0 ),
          ephemerisLineData_( )
    { }

    //! Default destructor.
    virtual ~ApproximatePlanetPositionsBase( ) { }

    //! Parse ephemeris line data.
    /*!
     * Parses ephemeris line data.
     * \param firstLineNumber First line number.
     */
    void parseEphemerisLineData_( const unsigned int& firstLineNumber );

    //! Parse line data for extra terms for ephemeris.
    /*!
     * Parses ephemeris line data for necessary extra terms.
     * \param lineNumber Line number.
     */
    void parseExtraTermsEphemerisLineData_( const unsigned int& lineNumber );

    //! Load in ephemeris data for planets.
    /*!
     * This method opens and parses the p_elem_t2.txt ephemeris files for the planet positions.
     * The resulting data is stored in
     * ApproximatePlanetPositionsBase::containerOfDataFromEphemerisFile_ to be used in the
     * generation of planet ephemeris. This method is automatically invoked if you call
     * ApproximatePlanetPositionsBase::setPlanet( BodiesWithEphemerisData ).
     */
    void reloadData( );

    //! Returns the gravitational parameter of the Sun.
    /*!
     *  Returns the gravitational parameter of the Sun that is used in the calculations.
     *  \return Gravitational parameter of the Sun.
     */
    double getSunGravitationalParameter( ){ return sunGravitationalParameter; }

protected:

    //! Set planet.
    /*!
     * Sets planet to retrieve ephemeris data for.
     * \param bodyWithEphemerisData Planet.
     */
    void setPlanet( BodiesWithEphemerisData bodyWithEphemerisData );

    //! Gravitational parameter of the Sun.
    /*!
     *  Gravitational parameter of the Sun.
     */
    const double sunGravitationalParameter;

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
     * by user when retrieving state.
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

    //! Keplerian elements of planet at given Julian date.
    /*!
     * Keplerian elements of planet at given Julian date.
     */
    Eigen::Vector6d planetKeplerianElementsAtGivenJulianDate_;

    //! String stream for ephemeris line data.
    /*!
     * String stream for ephemeris line data.
     */
    std::stringstream ephemerisLineData_;

private:
};

//! Typedef for shared-pointer to ApproximatePlanetPositionsBase object.
typedef std::shared_ptr< ApproximatePlanetPositionsBase > ApproximatePlanetPositionsBasePointer;

} // namespace ephemerides
} // namespace tudat

#endif // TUDAT_APPROXIMATE_PLANET_POSITIONS_BASE_H

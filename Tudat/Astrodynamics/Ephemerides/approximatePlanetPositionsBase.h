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
 *      110221    K. Kumar          File created.
 *      110224    K. Kumar          Renamed class and file.
 *      110629    L. van der Ham    Added function coplanar circular orbits.
 *      110803    L. van der Ham    Created base class and separated approximatePlanetPositions
 *                                  from approximatePlanetPositionsCircularCoplanar.
 *      110824    J. Leloux         Corrected doxygen documentation.
 *      120322    D. Dirkx          Modified to new Ephemeris interfaces.
 *      130121    K. Kumar          Added shared-ptr typedef.
 *      130120    D. Dirkx          Updated with new Julian day + seconds since Julian day input.
 *      140115    E. Brandon        Corrected doxygen documentation.
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the Major Planets,
 *          http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf, last accessed: 24 February, 2011.
 *
 *    Notes
 *
 */

#ifndef TUDAT_APPROXIMATE_PLANET_POSITIONS_BASE_H
#define TUDAT_APPROXIMATE_PLANET_POSITIONS_BASE_H

#include <cmath>
#include <fstream>
#include <map>
#include <string>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <TudatCore/Mathematics/BasicMathematics/basicMathematicsFunctions.h>

#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositionsDataContainer.h"
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

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
    basic_mathematics::Vector6d planetKeplerianElementsAtGivenJulianDate_;

    //! String stream for ephemeris line data.
    /*!
     * String stream for ephemeris line data.
     */
    std::stringstream ephemerisLineData_;

private:
};

//! Typedef for shared-pointer to ApproximatePlanetPositionsBase object.
typedef boost::shared_ptr< ApproximatePlanetPositionsBase > ApproximatePlanetPositionsBasePointer;

} // namespace ephemerides
} // namespace tudat

#endif // TUDAT_APPROXIMATE_PLANET_POSITIONS_BASE_H

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
 *      110224    K. Kumar          Creation of code.
 *      110810    J. Leloux         Corrected doxygen documentation.
 *      120322    D. Dirkx          Modified to new Ephemeris interfaces.
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the Major Planets,
 *          http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf, last accessed: 24 February, 2011.
 *
 */

#ifndef TUDAT_APPROXIMATE_PLANET_POSITIONS_DATA_CONTAINER_H
#define TUDAT_APPROXIMATE_PLANET_POSITIONS_DATA_CONTAINER_H

#include <iostream>
#include <string>

namespace tudat
{
namespace ephemerides
{

//! JPL "Approximate Positions of Major Planets" data container class.
/*!
 * Data container class for JPL "Approximate Positions of Major Planets" ephemeris data.
 */
struct ApproximatePlanetPositionsDataContainer
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    ApproximatePlanetPositionsDataContainer( )
        : semiMajorAxis_( -0.0 ),
          eccentricity_( -0.0 ),
          inclination_( -0.0 ),
          meanLongitude_( -0.0 ),
          longitudeOfPerihelion_( -0.0 ),
          longitudeOfAscendingNode_( -0.0 ),
          rateOfChangeOfSemiMajorAxis_( -0.0 ),
          rateOfChangeOfEccentricity_( -0.0 ),
          rateOfChangeOfInclination_( -0.0 ),
          rateOfChangeOfMeanLongitude_( -0.0 ),
          rateOfChangeOfLongitudeOfPerihelion_( -0.0 ),
          rateOfChangeOfLongitudeOfAscendingNode_( -0.0 ),
          additionalTermB_( -0.0 ),
          additionalTermC_( -0.0 ),
          additionalTermS_( -0.0 ),
          additionalTermF_( -0.0 )
    { }

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param approximatePlanetPositionsDataContainer Aproximate planet
     *          positions data container.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream,
                                     ApproximatePlanetPositionsDataContainer&
                                     approximatePlanetPositionsDataContainer )
    {
        using std::endl;

        stream << "This is an ApproximatePlanetPositionsDataContainer object. " << endl;
        stream << "The data corresponds to the table entry for "
               << approximatePlanetPositionsDataContainer.planetName_ << endl;
        stream << "The semi-major axis in AU is set to: "
               << approximatePlanetPositionsDataContainer.semiMajorAxis_ << endl;
        stream << "The eccentricity in radians is set to: "
               << approximatePlanetPositionsDataContainer.eccentricity_ << endl;
        stream << "The inclination in degrees is set to: "
               << approximatePlanetPositionsDataContainer.inclination_ << endl;
        stream << "The mean longitude in degrees is set to: "
               << approximatePlanetPositionsDataContainer.meanLongitude_ << endl;
        stream << "The longitude of perihelion in degrees is set to: "
               << approximatePlanetPositionsDataContainer.longitudeOfPerihelion_ << endl;
        stream << "The longitude of the ascending node in degrees is set to: "
               << approximatePlanetPositionsDataContainer.longitudeOfAscendingNode_<< endl;
        stream << "The rate of change of semi-major axis AU per century is set to: "
               << approximatePlanetPositionsDataContainer.rateOfChangeOfSemiMajorAxis_ << endl;
        stream << "The rate of change of eccentricity in radians per century is set to: "
               << approximatePlanetPositionsDataContainer.rateOfChangeOfEccentricity_ << endl;
        stream << "The rate of change of inclination in degrees per century is set to: "
               << approximatePlanetPositionsDataContainer.rateOfChangeOfInclination_ << endl;
        stream << "The rate of change of mean longitude in degrees per century is set to: "
               << approximatePlanetPositionsDataContainer.rateOfChangeOfMeanLongitude_ << endl;
        stream << "The rate of change of longitude of perihelion in degrees per century is set to: "
               << approximatePlanetPositionsDataContainer.rateOfChangeOfLongitudeOfPerihelion_ << endl;
        stream << "The rate of change of longitude of ascending node in degrees per century is set "
               << "to: " << approximatePlanetPositionsDataContainer
                  .rateOfChangeOfLongitudeOfAscendingNode_ << endl;

        // Check if additional terms are defined for outer planets.
        if ( approximatePlanetPositionsDataContainer.additionalTermB_ != -0.0 )
        {
            stream << "The additional term B is set to: "
                   << approximatePlanetPositionsDataContainer.additionalTermB_ << endl;
            stream << "The additional term C is set to: "
                   << approximatePlanetPositionsDataContainer.additionalTermC_ << endl;
            stream << "The additional term S is set to: "
                   << approximatePlanetPositionsDataContainer.additionalTermS_ << endl;
            stream << "The additional term F is set to: "
                   << approximatePlanetPositionsDataContainer .additionalTermF_ << endl;
        }

        // Return stream.
        return stream;
    }

    //! Planet name.
    /*!
     * Planet name
     */
    std::string planetName_;

    //! Semi-major axis.
    /*!
     * Semi-major axis given in Astronomical Units.
     */
    double semiMajorAxis_;

    //! Eccentricity.
    /*!
     * Eccentricity given in radians.
     */
    double eccentricity_;

    //! Inclination.
    /*!
     * Inclination given in radians.
     */
    double inclination_;

    //! Mean longitude.
    /*!
     * Mean longitude given in radians.
     */
    double meanLongitude_;

    //! Longitude of perihelion.
    /*!
     * Longitude of perihelion given in radians.
     */
    double longitudeOfPerihelion_;

    //! Longitude of ascending node.
    /*!
     * Longitude of the ascending node given in radians.
     */
    double longitudeOfAscendingNode_;

    //! Rate of change of semi-major axis.
    /*!
     * Rate of change of semi-major axis given in Astronomical Units per
     * century.
     */
    double rateOfChangeOfSemiMajorAxis_;

    //! Rate of change of eccentricity.
    /*!
     * Rate of change of eccentricity given in radians per century.
     */
    double rateOfChangeOfEccentricity_;

    //! Rate of change of inclination.
    /*!
     * Rate of change of inclination given in degrees per century.
     */
    double rateOfChangeOfInclination_;

    //! Rate of change of mean longitude.
    /*!
     * Rate of change of mean longitude given in degrees per century.
     */
    double rateOfChangeOfMeanLongitude_;

    //! Rate of change of longitude of perihelion.
    /*!
     * Rate of change of longitude of perihelion given in degrees per century.
     */
    double rateOfChangeOfLongitudeOfPerihelion_;

    //! Rate of change of longitude of ascending node.
    /*!
     * Rate of change of longitude of the ascending node given in degrees per
     * century.
     */
    double rateOfChangeOfLongitudeOfAscendingNode_;

    //! Additional term b.
    /*!
     * Additional term b, required for the computation of the mean anomaly for
     * Jupiter through Pluto, for the period 3000 BC to 3000 AD.
     */
    double additionalTermB_;

    //! Additional term c.
    /*!
     * Additional term c, required for the computation of the mean anomaly for
     * Jupiter through Pluto, for the period 3000 BC to 3000 AD.
     */
    double additionalTermC_;

    //! Additional term s.
    /*!
     * Additional term s, required for the computation of the mean anomaly for
     * Jupiter through Pluto, for the period 3000 BC to 3000 AD.
     */
    double additionalTermS_;

    //! Additional term f.
    /*!
     * Additional term f, required for the computation of the mean anomaly for
     * Jupiter through Pluto, for the period 3000 BC to 3000 AD.
     */
    double additionalTermF_;

protected:

private:
};

} // namespace ephemerides
} // namespace tudat

#endif // TUDAT_APPROXIMATE_PLANET_POSITIONS_DATA_CONTAINER_H

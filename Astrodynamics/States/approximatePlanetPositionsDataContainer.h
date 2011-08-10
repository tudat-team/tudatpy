/*! \file approximatePlanetPositionsDataContainer.h
 *    This header file contains the definition of a data container class for
 *    data extracted from the JPL "Approximate Positions of Major Planets"
 *    ephemeris for a specific planet ( http://ssd.jpl.nasa.gov/?planet_pos ).
 *    The ephemeris file used is for the period 3000 BC to 3000 AD.
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
 *    Date created      : 24 February, 2011
 *    Last modified     : 10 August, 2011
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
 *      110224    K. Kumar          First creation of code.
 *      110810    J. Leloux         Corrected doxygen documentation.
 */

#ifndef APPROXIMATEPLANETPOSITIONSDATACONTAINER_H
#define APPROXIMATEPLANETPOSITIONSDATACONTAINER_H

// Include statements.
#include <iostream>
#include <string>

// Using declarations.
using std::string;
using std::endl;

//! JPL "Approximate Positions of Major Planets" data container class.
/*!
 * Data container class for JPL "Approximate Positions of Major Planets"
 * ephemeris data.
 */
class ApproximatePlanetPositionsDataContainer
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    ApproximatePlanetPositionsDataContainer( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~ApproximatePlanetPositionsDataContainer( );

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
                                     approximatePlanetPositionsDataContainer );

    //! Planet name.
    /*!
     * Planet name
     */
    string planetName_;

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

#endif // APPROXIMATEPLANETPOSITIONSDATACONTAINER_H

// End of file.

/*! \file approximatePlanetPositions.h
 *    This header file contains the definition of an ephemeris class that makes
 *    use of the JPL "Approximate Positions of Major Planets"
 *    ( http://ssd.jpl.nasa.gov/?planet_pos ) to retrieve ephemeris data for a
 *    specific planet. The ephemeris file used is for the period 3000 BC to
 *    3000 AD.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 3
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
 *      110221    K. Kumar          First creation of code.
 *      110224    K. Kumar          Renamed class and file.
 *      110810    J. Leloux         Corrected doxygen documentation.
 */

#ifndef APPROXIMATEPLANETPOSITIONS_H
#define APPROXIMATEPLANETPOSITIONS_H

// Include statements.
#include "approximatePlanetPositionsBase.h"
#include "convertMeanAnomalyToEccentricAnomaly.h"
#include "newtonRaphson.h"

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
    ApproximatePlanetPositions( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~ApproximatePlanetPositions( );

    //! Get state from ephemeris.
    /*!
     * Returns state in Cartesian elements from ephemeris.
     * \return State in Cartesian elements from ephemeris.
     */
    CartesianElements* getStateFromEphemeris( const double& julianDate );

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

    //! Convert mean anomaly to eccentric anomaly.
    /*!
     * Convert mean anomaly to eccentric anomaly.
     */
    orbital_element_conversions::ConvertMeanAnomalyToEccentricAnomaly
            convertMeanAnomalyToEccentricAnomaly_;

    //! Newton-Raphson method.
    /*!
     * Newton-Raphson method.
     */
    NewtonRaphson newtonRaphson_;
};

#endif // APPROXIMATEPLANETPOSITIONS_H

// End of file.

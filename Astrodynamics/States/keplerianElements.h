/*! \file keplerianElements.h
 *    This header file contains the Keplerian elements class included in Tudat.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 20 October, 2010
 *    Last modified     : 2 December, 2010
 *
 *    References
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
 *      101020    K. Kumar          First creation of code.
 *      101022    K. Kumar          Added set/get functions and completed code
 *                                  comments.
 *      101110    K. Kumar          Added get functions for auxilliary
 *                                  parameters.
 *      101130    E. Iorfida        Added set function for semi-latus rectum.
 */

#ifndef KEPLERIANELEMENTS_H
#define KEPLERIANELEMENTS_H

// Include statements.
#include "orbitalElements.h"
#include "basicMathematicsFunctions.h"

//! Keplerian elements class.
/*!
 * Keplerian elements class.
 */
class KeplerianElements : public OrbitalElements
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    KeplerianElements( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~KeplerianElements( );

    //! Set semi-major axis.
    /*!
     * Sets semi-major axis.
     * \param semiMajorAxis Semi-major axis.
     */
    void setSemiMajorAxis( const double& semiMajorAxis );

    //! Set eccentricity.
    /*!
     * Sets eccentricity.
     * \param eccentricity Eccentricity.
     */
    void setEccentricity( const double& eccentricity );

    //! Set inclination.
    /*!
     * Sets inclination.
     * \param inclination Inclination.
     */
    void setInclination( const double& inclination );

    //! Set argument of periapsis.
    /*!
     * Sets argument of periapsis.
     * \param argumentOfPeriapsis Argument of periapsis.
     */
    void setArgumentOfPeriapsis( const double& argumentOfPeriapsis );

    //! Sets right ascension of ascending node.
    /*!
     * Set right ascension of ascending node.
     * \param rightAscensionOfAscendingNode Right ascension of ascending node.
     */
    void setRightAscensionOfAscendingNode( const double&
                                           rightAscensionOfAscendingNode );

    //! Set true anomaly.
    /*!
     * Sets true anomaly.
     * \param trueAnomaly True anomaly.
     */
    void setTrueAnomaly( const double& trueAnomaly );

    //! Set semi-latus rectum ( for parabolic orbits ).
    /*!
     * This function sets the semi-latus rectum for parabolic orbits.
     * The semi-latus rectum is not computed from the orbital parameters
     * since for a parabola the semi-major axis is undefined. This function
     * must only be used in conjunction with parabolic orbits.
     * \param semiLatusRectum Semi-latus rectum.
     */
    void setSemiLatusRectum( const double& semiLatusRectum );

    //! Get semi-major axis.
    /*!
     * Get semi-major axis.
     * \return Semi-major axis.
     */
    double& getSemiMajorAxis( );

    //! Get eccentricity.
    /*!
     * Get eccentricity.
     * \return Eccentricity.
     */
    double& getEccentricity( );

    //! Get inclination.
    /*!
     * Get inclination.
     * \return Inclination.
     */
    double& getInclination( );

    //! Get argument of periapsis.
    /*!
     * Get argument of periapsis.
     * \return Argument of periapsis.
     */
    double& getArgumentOfPeriapsis( );

    //! Get right ascension of ascending node.
    /*!
     * Get right ascension of ascending node.
     * \return Right ascension of ascending node.
     */
    double& getRightAscensionOfAscendingNode( );

    //! Get true anomaly.
    /*!
     * Get true anomaly.
     * \return True anomaly.
     */
    double& getTrueAnomaly( );

    //! Get longitude of periapsis.
    /*!
     * Get longitude of periapsis.
     * \return Longitude of periapsis.
     */
    double getLongitudeOfPeriapsis( );

    //! Get true longitude.
    /*!
     * Get true longitude.
     * \return True longitude.
     */
    double getTrueLongitude( );

    //! Get argument of latitude.
    /*!
     * Get argument of latitude.
     * \return Argument of latitude.
     */
    double getArgumentOfLatitude( );

    //! Get semi-latus rectum.
    /*!
     * Get semi-latus rectum.
     * \return Semi-latus rectum.
     */
    double getSemiLatusRectum( );

protected:

private:

    //! Semi-latus rectum ( for parabolas ).
    /*!
     * Semi-latus rectum, only applicable for parabolas.
     */
    double semiLatusRectum_;
};

#endif // KEPLERIANELEMENTS_H

// End of file.

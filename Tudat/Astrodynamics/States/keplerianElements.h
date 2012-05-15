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
 *      101020    K. Kumar          First creation of code.
 *      101022    K. Kumar          Added set/get functions and completed code comments.
 *      101110    K. Kumar          Added get functions for auxilliary parameters.
 *      101130    E. Iorfida        Added set function for semi-latus rectum.
 *      110310    K. Kumar          Changed right ascension of ascending node to longitude of
 *                                  ascending node.
 *      120511    K. Kumar          Added enum for Keplerian element indices.
 *
 *    References
 *
 */

#ifndef TUDAT_KEPLERIAN_ELEMENTS_H
#define TUDAT_KEPLERIAN_ELEMENTS_H

#include "Tudat/Astrodynamics/States/state.h"

#include <TudatCore/Mathematics/BasicMathematics/basicMathematicsFunctions.h>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

namespace tudat
{
namespace astrodynamics
{
namespace states
{

//! Keplerian element indices.
/*!
 * Keplerian element vector indices.
 */
enum KeplerianElementsIndices
{
    semiMajorAxisIndex,
    eccentricityIndex,
    inclinationIndex,
    argumentOfPeriapsisIndex,
    longitudeOfAscendingNodeIndex,
    trueAnomalyIndex,
    semiLatusRectumIndex = 0,
    longitudeOfPeriapsisIndex = 1,
    trueLongitudeIndex = 2,
    argumentOfLatitudeIndex = 3
};

//! Keplerian elements class.
/*!
 * Keplerian elements class.
 */
class KeplerianElements : public State
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    KeplerianElements( ) { state.setZero( 6 ); }

    //! Set semi-major axis.
    /*!
     * Sets semi-major axis.
     * \param semiMajorAxis Semi-major axis.
     */
    void setSemiMajorAxis( const double semiMajorAxis )
    {
        state( semiMajorAxisIndex ) = semiMajorAxis;
    }

    //! Set eccentricity.
    /*!
     * Sets eccentricity.
     * \param eccentricity Eccentricity.
     */
    void setEccentricity( const double eccentricity )
    {
        state( eccentricityIndex ) = eccentricity;
    }

    //! Set inclination.
    /*!
     * Sets inclination.
     * \param inclination Inclination.
     */
    void setInclination( const double inclination ) { state( inclinationIndex ) = inclination; }

    //! Set argument of periapsis.
    /*!
     * Sets argument of periapsis.
     * \param argumentOfPeriapsis Argument of periapsis.
     */
    void setArgumentOfPeriapsis( const double argumentOfPeriapsis )
    {
        state( argumentOfPeriapsisIndex ) = argumentOfPeriapsis;
    }

    //! Sets longitude of ascending node.
    /*!
     * Set longitude of ascending node.
     * \param longitudeOfAscendingNode Longitude of ascending node.
     */
    void setLongitudeOfAscendingNode( const double longitudeOfAscendingNode )
    {
        state( longitudeOfAscendingNodeIndex ) = longitudeOfAscendingNode;
    }

    //! Set true anomaly.
    /*!
     * Sets true anomaly.
     * \param trueAnomaly True anomaly.
     */
    void setTrueAnomaly( const double trueAnomaly ) { state( trueAnomalyIndex ) = trueAnomaly; }

    //! Set semi-latus rectum ( for parabolic orbits ).
    /*!
     * Sets the semi-latus rectum for parabolic orbits. The semi-latus rectum is not
     * computed from the orbital parameters since for a parabola the semi-major axis is undefined.
     * This function must only be used in conjunction with parabolic orbits.
     * \param semiLatusRectum Semi-latus rectum.
     */
    void setSemiLatusRectum( const double semiLatusRectum )
    {
        state( semiLatusRectumIndex ) = semiLatusRectum;
    }

    //! Get semi-major axis.
    /*!
     * Get semi-major axis.
     * \return Semi-major axis.
     */
    double getSemiMajorAxis( ) { return state( semiMajorAxisIndex ); }

    //! Get eccentricity.
    /*!
     * Get eccentricity.
     * \return Eccentricity.
     */
    double getEccentricity( ) { return state( eccentricityIndex ); }

    //! Get inclination.
    /*!
     * Get inclination.
     * \return Inclination.
     */
    double getInclination( ) { return state( inclinationIndex ); }

    //! Get argument of periapsis.
    /*!
     * Get argument of periapsis.
     * \return Argument of periapsis.
     */
    double getArgumentOfPeriapsis( ) { return state( argumentOfPeriapsisIndex ); }

    //! Get longitude of ascending node.
    /*!
     * Get longitude of ascending node.
     * \return Longitude of ascending node.
     */
    double getLongitudeOfAscendingNode( ) { return state( longitudeOfAscendingNodeIndex ); }

    //! Get true anomaly.
    /*!
     * Get true anomaly.
     * \return True anomaly.
     */
    double getTrueAnomaly( ) { return state( trueAnomalyIndex ); }

    //! Get longitude of periapsis.
    /*!
     * Get longitude of periapsis.
     * \return Longitude of periapsis.
     */
    double getLongitudeOfPeriapsis( )
    {
        return getArgumentOfPeriapsis( ) + getLongitudeOfAscendingNode( );
    }

    //! Get true longitude.
    /*!
     * Get true longitude.
     * \return True longitude.
     */
    double getTrueLongitude( )
    {
        return getArgumentOfPeriapsis( ) + getLongitudeOfAscendingNode( ) + getTrueAnomaly( );
    }

    //! Get argument of latitude.
    /*!
     * Get argument of latitude.
     * \return Argument of latitude.
     */
    double getArgumentOfLatitude( )
    {
        return getArgumentOfPeriapsis( ) + getTrueAnomaly( );
    }

    //! Get semi-latus rectum.
    /*!
     * Get semi-latus rectum This only works if it has been set, which is necessary for a parabola.
     * It does not compute the semi-latus rectum from the semi-major axis and the eccentricity,
     * since in the case of a parabola, the semi-major axis is not defined.
     * \return Semi-latus rectum.
     */
    double getSemiLatusRectum( ) { return state( semiLatusRectumIndex ); }

protected:

private:
};

} // namespace states
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_KEPLERIAN_ELEMENTS_H

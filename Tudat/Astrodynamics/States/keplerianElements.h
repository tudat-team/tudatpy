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
 *      101020    K. Kumar          Creation of code.
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

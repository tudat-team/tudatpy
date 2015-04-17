/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      101020    K. Kumar          First creation of code.
 *      101025    E. Iorfida        Fulfillment of the code with gravitational parameter.
 *      101130    E. Iorfida        Gravitational parameter removed.
 *      101202    J. Melman         Replaced endif statement and changed. Doxygen return statement.
 *      101203    E. Iorfida        Gravitational parameter added.
 *      101219    J. Melman         Doxygen comment on gravitational parameter added.
 *      110128    K. Kumar          Changed references to pointers for functions.
 *      110510    K. Kumar          Updated conversion functions to not use dynamic memory
 *                                  allocation.
 *      110805    K. Kumar          Added mean motion to semi-major axis conversion.
 *      110810    J. Leloux         Corrected doxygen documentation.
 *      120131    K. Kumar          Adapted for Tudat Core, interfaces changed to use VectorXd,
 *                                  only Keplerian <-> Cartesian conversions included.
 *      120206    K. Kumar          Added wrapper functions for orbital element conversions when
 *                                  eccentricity is not known a priori (if-statement to choose
 *                                  between elliptical and hyperbolic orbits).
 *      120422    K. Kumar          Added Doxygen notes for Cartesian -> Keplerian conversion.
 *      121205    K. Kumar          Migrated namespace to directory-based protocol and added
 *                                  backwards compatibility.
 *
 *    References
 *      Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series, VA, 2002.
 *      Wertz, J. R. Mission geometry; orbit and constellation design and management.
 *      Mengali, G., Quarta, A.A. Fondamenti di meccanica del volo spaziale.
 *      Wertz, J.R. Mission Geometry; Orbit and Constellation Design and Management, Spacecraft
 *          Orbit and Attitude Systems, Microcosm Press, Kluwer Academic Publishers, 2001.
 *      Advanced Concepts Team, ESA. Keplerian Toolbox, http://sourceforge.net/projects/keptoolbox,
 *          last accessed: 21st April, 2012.
 *
 *    Notes
 *      Backwards compatibility of namespaces is implemented for Tudat Core 2 in this file. The
 *      code block marked "DEPRECATED!" at the end of the file should be removed in Tudat Core 3.
 *
 */

#ifndef TUDAT_ORBITAL_ELEMENT_CONVERSIONS_H
#define TUDAT_ORBITAL_ELEMENT_CONVERSIONS_H

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"

namespace tudat
{

namespace orbital_element_conversions
{

//! Convert Keplerian to Cartesian orbital elements.
/*!
 * Converts Keplerian to Cartesian orbital elements (Chobotov, 2002). Use the 
 * CartesianElementVectorIndices enum to access the individual orbital element components in the 
 * storage vector.
 *
 * \param keplerianElements Vector containing Keplerian elements.                         \n
 *          <em>  
 *                          Order of elements is important! \n            
 *                          keplerianElements( 0 ) = semiMajorAxis,                   [m] \n
 *                          keplerianElements( 1 ) = eccentricity,                    [-] \n
 *                          keplerianElements( 2 ) = inclination,                   [rad] \n
 *                          keplerianElements( 3 ) = argument of periapsis,         [rad] \n
 *                          keplerianElements( 4 ) = longitude of ascending node,   [rad] \n
 *                          keplerianElements( 5 ) = true anomaly.                  [rad] \n 
 *          </em>
 *        WARNING: If eccentricity is 1.0 within machine precision, 
 *        keplerianElements( 0 ) = semi-latus rectum.      
 *                                    
 * \param centralBodyGravitationalParameter Gravitational parameter of central body [m^3 s^-2]. 
 *   
 * \return Converted state in Cartesian elements.                         \n
 *         <em>  
 *         Order of elements is important!                                \n
 *         cartesianElements( 0 ) = x-position coordinate,            [m] \n
 *         cartesianElements( 1 ) = y-position coordinate,            [m] \n
 *         cartesianElements( 2 ) = z-position coordinate,            [m] \n
 *         cartesianElements( 3 ) = x-velocity coordinate,          [m/s] \n
 *         cartesianElements( 4 ) = y-velocity coordinate,          [m/s] \n
 *         cartesianElements( 5 ) = z-velocity coordinate.          [m/s] \n 
 *         </em>
 *
 * \sa CartesianElementVectorIndices()
 */
basic_mathematics::Vector6d convertKeplerianToCartesianElements(
        const basic_mathematics::Vector6d& keplerianElements,
        const double centralBodyGravitationalParameter );

//! Convert Cartesian to Keplerian orbital elements.
/*!
 * Converts Cartesian to Keplerian orbital elements.
 * \param cartesianElements Vector containing Cartesian elements. Order of elements is important!
 *          cartesianElements( 0 ) = x-position coordinate,                                     [m]
 *          cartesianElements( 1 ) = y-position coordinate,                                     [m]
 *          cartesianElements( 2 ) = z-position coordinate,                                     [m]
 *          cartesianElements( 3 ) = x-velocity coordinate,                                   [m/s]
 *          cartesianElements( 4 ) = y-velocity coordinate,                                   [m/s]
 *          cartesianElements( 5 ) = z-velocity coordinate.                                   [m/s]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.
 * \return Converted state in Keplerian elements. The order of elements is fixed!
 *          keplerianElements( 0 ) = semiMajorAxis,                                             [m]
 *          keplerianElements( 1 ) = eccentricity,                                              [-]
 *          keplerianElements( 2 ) = inclination,                                             [rad]
 *          keplerianElements( 3 ) = argument of periapsis,                                   [rad]
 *          keplerianElements( 4 ) = longitude of ascending node,                             [rad]
 *          keplerianElements( 5 ) = true anomaly.                                            [rad]
 *          WARNING: If eccentricity is 1.0 within 1.0e-15,
 *          keplerianElements( 0 ) = semi-latus rectum, since the orbit is parabolic.
 *          WARNING: If eccentricity is 0.0 within 1.0e-15,
 *          argument of periapsis is set to 0.0, since the orbit is circular.
 *          WARNING: If inclination is 0.0 within 1.0e-15,
 *          longitude of ascending node is set to 0.0, since the orbit is equatorial.
 *          The tolerance 1.0e-15 is hard-coded, as it should not be changed for performance
 *          reasons, unless required for specific scenarios. In those cases, users are expected
 *          to update the internal tolerance to the required value. Below these tolerance values
 *          for eccentricity and inclination, the orbit is considered to be a limit case.
 *          Essentially, special solutions are then used for parabolic, circular inclined,
 *          non-circular equatorial, and circular equatorial orbits. These special solutions are
 *          required because of singularities in the classical Keplerian elements. If high
 *          precision is required near these singularities, users are encouraged to consider using
 *          other elements, such as modified equinoctial elements. It should be noted that
 *          modified equinoctial elements also suffer from singularities, but not for zero
 *          eccentricity and inclination.
 */
basic_mathematics::Vector6d convertCartesianToKeplerianElements(
        const basic_mathematics::Vector6d& cartesianElements,
        const double centralBodyGravitationalParameter );

//! Convert true anomaly to (elliptical) eccentric anomaly.
/*!
 * Converts true anomaly to eccentric anomaly for elliptical orbits ( 0 <= eccentricity < 1.0 ).
 * The equations used can be found in (Chobotov, 2002).
 * \param trueAnomaly True anomaly.                                                           [rad]
 * \param eccentricity Eccentricity.                                                            [-]
 * \return (Elliptical) Eccentric anomaly.                                                    [rad]
 */
double convertTrueAnomalyToEllipticalEccentricAnomaly( const double trueAnomaly,
                                                       const double eccentricity );

//! Convert true anomaly to hyperbolic eccentric anomaly.
/*!
 * Converts true anomaly to hyperbolic eccentric anomaly for hyperbolic orbits
 * ( eccentricity > 1.0 ). The equations used can be found in (Chobotov, 2002).
 * \param trueAnomaly True anomaly.                                                           [rad]
 * \param eccentricity Eccentricity.                                                            [-]
 * \return Hyperbolic eccentric anomaly.                                                      [rad]
 */
double convertTrueAnomalyToHyperbolicEccentricAnomaly( const double trueAnomaly,
                                                       const double eccentricity );

//! Convert true anomaly to eccentric anomaly.
/*!
 * Converts true anomaly to eccentric anomaly for elliptical and hyperbolic orbits
 * ( eccentricity < 1.0 && eccentricity > 1.0 ). This function is essentially a wrapper for
 * convertTrueAnomalyToEllipticalEccentricAnomaly() and
 * convertTrueAnomalyToHyperbolicEccentricAnomaly(). It should be used in cases where the
 * eccentricity of the orbit is not known a priori. Currently, this implementation performs a
 * check on the eccentricity and throws an error for eccentricity < 0.0 and parabolic orbits, which
 * have not been implemented. The equations used can be found in (Chobotov, 2002).
 * \param eccentricAnomaly Eccentric anomaly.                                                 [rad]
 * \param eccentricity Eccentricity.                                                            [-]
 * \return True anomaly.                                                                      [rad]
 */
double convertTrueAnomalyToEccentricAnomaly( const double eccentricAnomaly,
                                             const double eccentricity );

//! Convert (elliptical) eccentric anomaly to true anomaly.
/*!
 * Converts eccentric anomaly to true anomaly for elliptical orbits ( 0 <= eccentricity < 1.0 ).
 * The equations used can be found in (Chobotov, 2002).
 * \param ellipticalEccentricAnomaly Elliptical eccentric anomaly.                            [rad]
 * \param eccentricity Eccentricity.                                                            [-]
 * \return True anomaly.                                                                      [rad]
 */
double convertEllipticalEccentricAnomalyToTrueAnomaly( const double ellipticalEccentricAnomaly,
                                                       const double eccentricity );

//! Convert hyperbolic eccentric anomaly to true anomaly.
/*!
 * Converts hyperbolic eccentric anomaly to true anomaly for hyperbolic orbits
 * ( eccentricity > 1.0 ). The equations used can be found in (Chobotov, 2002).
 * \param hyperbolicEccentricAnomaly Hyperbolic eccentric anomaly.                            [rad]
 * \param eccentricity Eccentricity.                                                            [-]
 * \return True anomaly.                                                                      [rad]
 */
double convertHyperbolicEccentricAnomalyToTrueAnomaly( const double hyperbolicEccentricAnomaly,
                                                       const double eccentricity );

//! Convert eccentric anomaly to true anomaly.
/*!
 * Converts eccentric anomaly to true anomaly for elliptical and hyperbolic orbits
 * ( eccentricity < 1.0 && eccentricity > 1.0 ). This function is essentially a wrapper for
 * convertEllipticalEccentricAnomalyToTrueAnomaly() and
 * convertHyperbolicEccentricAnomalyToTrueAnomaly(). It should be used in cases where the
 * eccentricity of the orbit is not known a priori. Currently, this implementation performs a
 * check on the eccentricity and throws an error for eccentricity < 0.0 and parabolic orbits, which
 * have not been implemented. The equations used can be found in (Chobotov, 2002).
 * \param eccentricAnomaly Eccentric anomaly.                                                 [rad]
 * \param eccentricity Eccentricity.                                                            [-]
 * \return True anomaly.                                                                      [rad]
 */
double convertEccentricAnomalyToTrueAnomaly( const double eccentricAnomaly,
                                             const double eccentricity );

//! Convert (elliptical) eccentric anomaly to mean anomaly.
/*!
 * Converts eccentric anomaly to mean anomaly for elliptical orbits ( 0 <= eccentricity < 1.0 ).
 * The equations used can be found in (Chobotov, 2002).
 * \param eccentricity Eccentricity.                                                            [-]
 * \param ellipticalEccentricAnomaly (Elliptical) eccentric anomaly [rad].
 * \return Mean anomaly [rad].
 */
double convertEllipticalEccentricAnomalyToMeanAnomaly( const double ellipticalEccentricAnomaly,
                                                       const double eccentricity );

//! Convert hyperbolic eccentric anomaly to mean anomaly.
/*!
 * Converts hyperbolic eccentric anomaly to mean anomaly for hyperbolic orbits
 * ( eccentricity > 1.0 ). The equations used can be found in (Chobotov, 2002).
 * \param hyperbolicEccentricAnomaly Hyperbolic eccentric anomaly.                            [rad]
 * \param eccentricity Eccentricity.                                                            [-]
 * \return Mean anomaly.                                                                      [rad]
 */
double convertHyperbolicEccentricAnomalyToMeanAnomaly( const double hyperbolicEccentricAnomaly,
                                                       const double eccentricity );

//! Convert eccentric anomaly to mean anomaly.
/*!
 * Converts eccentric anomaly to mean anomaly for elliptical and hyperbolic orbits
 * ( eccentricity < 1.0 && eccentricity > 1.0 ). This function is essentially a wrapper for
 * convertEllipticalEccentricAnomalyToMeanAnomaly() and
 * convertHyperbolicEccentricAnomalyToMeanAnomaly(). It should be used in cases where the
 * eccentricity of the orbit is not known a priori. Currently, this implementation performs a
 * check on the eccentricity and throws an error for eccentricity < 0.0 and parabolic orbits, which
 * have not been implemented. The equations used can be found in (Chobotov, 2002).
 * \param eccentricity Eccentricity.                                                            [-]
 * \param eccentricAnomaly Eccentric anomaly.                                                 [rad]
 * \return Mean anomaly.                                                                      [rad]
 */
double convertEccentricAnomalyToMeanAnomaly( const double eccentricAnomaly,
                                             const double eccentricity );

//! Convert elapsed time to (elliptical) mean anomaly change.
/*!
 * Converts elapsed time to mean anomaly change for elliptical orbits ( 0 <= eccentricity < 1.0 ).
 * The semi-major axis must be non-negative; this function will throw an error to indicate if the
 * semi-major axis is invalid. The equation used can be found in (Chobotov, 2002).
 * \param elapsedTime Elapsed time.                                                             [s]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \param semiMajorAxis Semi-major axis.                                                        [m]
 * \return (Elliptical) Mean anomaly change.                                                  [rad]
 */
double convertElapsedTimeToEllipticalMeanAnomalyChange(
        const double elapsedTime, const double centralBodyGravitationalParameter,
        const double semiMajorAxis );

//! Convert elapsed time to mean anomaly change for hyperbolic orbits.
/*!
 * Converts elapsed time to mean anomaly change for hyperbolic orbits ( eccentricity > 1.0 ).
 * The semi-major axis must be non-positive; this function will throw an error to indicate if the
 * semi-major axis is invalid. The equation used can be found in (Chobotov, 2002).
 * \param elapsedTime Elapsed time.                                                             [s]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \param semiMajorAxis Semi-major axis.                                                        [m]
 * \return Mean anomaly change.                                                               [rad]
 */
double convertElapsedTimeToHyperbolicMeanAnomalyChange(
        const double elapsedTime, const double centralBodyGravitationalParameter,
        const double semiMajorAxis );

//! Convert elapsed time to mean anomaly change.
/*!
 * Converts elapsed time to mean anomaly change for elliptical and hyperbolic orbits
 * ( eccentricity < 1.0 && eccentricity > 1.0 ). This function is essentially a wrapper for
 * convertElapsedTimeToEllipticalMeanAnomalyChange() and
 * convertElapsedTimeToHyperbolicMeanAnomalyChange(). It should be used in cases where the
 * eccentricity of the orbit is not known a priori. The equations used can be found in
 * (Wertz, 2001).
 * \param elapsedTime Elapsed time.                                                             [s]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \param semiMajorAxis Semi-major axis.                                                        [m]
 * \return Mean anomaly change.                                                               [rad]
 */
double convertElapsedTimeToMeanAnomalyChange(
        const double elapsedTime, const double centralBodyGravitationalParameter,
        const double semiMajorAxis );

//! Convert (elliptical) mean anomaly change to elapsed time.
/*!
 * Converts mean anomaly change to elapsed time for elliptical orbits ( 0 <= eccentricity < 1.0 ).
 * The equation used can be found in (Wertz, 2001). This function checks if the semi-major axis is
 * non-negative and throws an error if it not.
 * \param ellipticalMeanAnomalyChange (Elliptical) Mean anomaly change.                       [rad]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \param semiMajorAxis Semi-major axis.                                                        [m]
 * \return Elapsed time.                                                                        [s]
 */
double convertEllipticalMeanAnomalyChangeToElapsedTime(
        const double ellipticalMeanAnomalyChange, const double centralBodyGravitationalParameter,
        const double semiMajorAxis );

//! Convert hyperbolic mean anomaly change to elapsed time.
/*!
 * Converts mean anomaly change to elapsed time for hyperbolic orbits ( eccentricity > 1.0 ).
 * The equation used can be found in (Wertz, 2001). This function checks if the semi-major axis is
 * non-positive and throws an error if it not.
 * \param hyperbolicMeanAnomalyChange Hyperbolic mean anomaly change.                         [rad]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \param semiMajorAxis Semi-major axis.                                                        [m]
 * \return Elapsed time.                                                                        [s]
 */
double convertHyperbolicMeanAnomalyChangeToElapsedTime(
        const double hyperbolicMeanAnomalyChange, const double centralBodyGravitationalParameter,
        const double semiMajorAxis );

//! Convert mean anomaly change to elapsed time.
/*!
 * Converts mean anomaly change to elapsed time for elliptical and hyperbolic orbits
 * ( eccentricity < 1.0 && eccentricity > 1.0 ). This function is essentially a wrapper for
 * convertEllipticalMeanAnomalyChangeToElapsedTime() and
 * convertHyperbolicMeanAnomalyChangeToElapsedTime(). It should be used in cases where the
 * eccentricity of the orbit is not known a priori. The equations used can be found in
 * (Wertz, 2001).
 * \param meanAnomalyChange Mean anomaly change.                                              [rad]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \param semiMajorAxis Semi-major axis.                                                        [m]
 * \return Elapsed time.                                                                        [s]
 */
double convertMeanAnomalyChangeToElapsedTime(
        const double meanAnomalyChange, const double centralBodyGravitationalParameter,
        const double semiMajorAxis );

//! Convert (elliptical) mean motion to semi-major axis.
/*!
 * Converts mean motion to semi-major axis for elliptical orbits.
 * \param ellipticalMeanMotion (Elliptical) Mean motion.                                    [rad/s]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \return semiMajorAxis Semi-major axis.                                                       [m]
 */
double convertEllipticalMeanMotionToSemiMajorAxis(
        const double ellipticalMeanMotion, const double centralBodyGravitationalParameter );

//! Convert semi-major axis to elliptical mean motion.
/*!
 * Converts semi-major axis to elliptical mean motion.
 * \param semiMajorAxis Semi-major axis.                                                        [m]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \return ellipticalMeanMotion (Elliptical) Mean motion.                                   [rad/s]
 */
double convertSemiMajorAxisToEllipticalMeanMotion(
        const double semiMajorAxis, const double centralBodyGravitationalParameter );

} // namespace orbital_element_conversions

} // namespace tudat

#endif // TUDAT_ORBITAL_ELEMENT_CONVERSIONS_H

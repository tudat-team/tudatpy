/*! \file orbitalElementConversions.h
 *    This header file contains a namespace with orbital element conversion
 *    functions.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 10
 *    Check status      : Checked
 *
 *    Author            : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 20 October, 2010
 *    Last modified     : 10 August, 2011
 *
 *    References
 *      Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series,
 *          VA, 2002.
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      101025    E. Iorfida        Fulfillment of the code with gravitational
 *                                  parameter.
 *      101130    E. Iorfida        Gravitational parameter removed.
 *      101202    J. Melman         Replaced endif statement and changed.
 *                                  Doxygen return statement.
 *      101203    E. Iorfida        Gravitational parameter added.
 *      101219    J. Melman         Doxygen comment on gravitational parameter
 *                                  added.
 *      110128    K. Kumar          Changed references to pointers for
 *                                  functions.
 *      110510    K. Kumar          Updated conversion functions to not use
 *                                  dynamic memory allocation.
 *      110805    K. Kumar          Added mean motion to semi-major axis
 *                                  conversion.
 *      110810    J. Leloux         Corrected doxygen documentation.
 */

// Include statements.
#include "Astrodynamics/Bodies/CelestialBodies/celestialBody.h"
#include "Astrodynamics/States/cartesianElements.h"
#include "Astrodynamics/States/keplerianElements.h"
#include "Mathematics/LinearAlgebra/linearAlgebra.h"

#ifndef ORBITALELEMENTCONVERSIONS_H
#define ORBITALELEMENTCONVERSIONS_H

//! Orbital element conversions namespace.
/*!
 *  Orbital element conversions namespace.
 */
namespace orbital_element_conversions
{

//! Convert Keplerian to Cartesian orbital elements.
/*!
 * Converts Keplerian to Cartesian orbital elements.
 * \param pointerToKeplerianElements Pointer to KeplerianElements object.
 * \param pointerToCelestialBody Pointer to CelestialBody object.
 * \return CartesianElements object.
 */
CartesianElements convertKeplerianToCartesianElements(
        KeplerianElements* pointerToKeplerianElements,
        CelestialBody* pointerToCelestialBody );

//! Convert Cartesian to Keplerian orbital elements.
/*!
 * Converts Cartesian to Keplerian orbital elements.
 * \param pointerToCartesianElements Pointer to CartesianElements object.
 * \param pointerToCelestialBody Pointer to CelestialBody object.
 * \return KeplerianElements object.
 */
KeplerianElements convertCartesianToKeplerianElements(
        CartesianElements* pointerToCartesianElements,
        CelestialBody* pointerToCelestialBody );

//! Convert true anomaly to eccentric anomaly.
/*!
 * Converts true anomaly to eccentric anomaly for elliptical orbits
 * ( 0 <= eccentricity < 1.0 ). The equations used can be found in
 * ( Chobotov, 2000 ).
 * \param trueAnomaly True anomaly.
 * \param eccentricity Eccentricity.
 * \return Eccentric anomaly.
 */
double convertTrueAnomalyToEccentricAnomaly( const double& trueAnomaly,
                                             const double& eccentricity );

//! Convert eccentric anomaly to true anomaly.
/*!
 * Converts eccentric anomaly to true anomaly for elliptical orbits
 * ( 0 <= eccentricity < 1.0 ). The equations used can be found in
 * ( Chobotov, 2000 ).
 * \param eccentricAnomaly Eccentric anomaly.
 * \param eccentricity Eccentricity.
 * \return True anomaly.
 */
double convertEccentricAnomalyToTrueAnomaly( const double& eccentricAnomaly,
                                             const double& eccentricity );

//! Convert true anomaly to hyperbolic eccentric anomaly.
/*!
 * Converts true anomaly to hyperbolic eccentric anomaly for hyperbolic orbits
 * ( eccentricity > 1.0 ). The equations used can be found in
 * ( Chobotov, 2000 ).
 * \param trueAnomaly True anomaly.
 * \param eccentricity Eccentricity.
 * \return Hyperbolic eccentric anomaly.
 */
double convertTrueAnomalyToHyperbolicEccentricAnomaly(
        const double& trueAnomaly, const double& eccentricity );

//! Convert hyperbolic eccentric anomaly to true anomaly.
/*!
 * Converts hyperbolic eccentric anomaly to true anomaly for hyperbolic orbits
 * ( eccentricity > 1.0 ). The equations used can be found in
 * ( Chobotov, 2000 ).
 * \param hyperbolicEccentricAnomaly Hyperbolic eccentric anomaly.
 * \param eccentricity Eccentricity.
 * \return True anomaly.
 */
double convertHyperbolicEccentricAnomalyToTrueAnomaly(
        const double& hyperbolicEccentricAnomaly, const double& eccentricity );

//! Convert eccentric anomaly to mean anomaly.
/*!
 * Converts eccentric anomaly to mean anomaly for elliptical orbits
 * ( 0 <= eccentricity < 1.0 ). The equations used can be found in
 * ( Chobotov, 2000 ).
 * \param eccentricity Eccentricity.
 * \param eccentricAnomaly Eccentric anomaly.
 * \return Mean anomaly.
 */
double convertEccentricAnomalyToMeanAnomaly( const double& eccentricAnomaly,
                                             const double& eccentricity );

//! Class to convert mean anomaly to eccentric anomaly.
/*!
 * This class is defines and implements conversion from mean anomaly to
 * eccentric anomaly for elliptical orbits. The conversion does not work
 * for near-parabolic orbits at present, hence the range of eccentricity is
 * restricted to: 0 <= eccentricity < 0.8.
 */
class convertMeanAnomalyToEccentricAnomaly;

//! Convert hyperbolic eccentric anomaly to mean anomaly.
/*!
 * Converts hyperbolic eccentric anomaly to mean anomaly for hyperbolic orbits
 * ( eccentricity > 1.0 ). The equations used can be found in
 * ( Chobotov, 2000 ).
 * \param hyperbolicEccentricAnomaly Hyperbolic eccentric anomaly.
 * \param eccentricity Eccentricity.
 * \return Mean anomaly.
 */
double convertHyperbolicEccentricAnomalyToMeanAnomaly(
        const double& hyperbolicEccentricAnomaly, const double& eccentricity );

//! Class to convert mean anomaly to hyperbolic eccentric anomaly.
/*!
 * This class is defines and implements conversion from mean anomaly to
 * hyperbolic eccentric anomaly for hyperbolic orbits. The conversion does
 * not work for near-parabolic orbits at present, hence the range of
 * eccentricity is restricted to: eccentricity > 1.2.
 */
class convertMeanAnomalyToHyperbolicEccentricAnomaly;

//! Convert elapsed time to mean anomaly for elliptical orbits.
/*!
 * Converts elapsed time to mean anomaly for elliptical orbits
 * ( 0 <= eccentricity < 1.0 ). The equation used can be found in
 * ( Chobotov, 2000 ).
 * \param elapsedTime Elapsed time.
 * \param pointerToCentralBody Pointer to central body.
 * \param semiMajorAxis Semi-major axis.
 * \return Mean anomaly.
 */
double convertElapsedTimeToMeanAnomalyForEllipticalOrbits(
        const double& elapsedTime, CelestialBody* pointerToCentralBody,
        const double& semiMajorAxis );

//! Convert mean anomaly to elapsed time.
/*!
 * Converts mean anomaly to elapsed time for elliptical orbits
 * ( 0 <= eccentricity < 1.0 ). The equation used can be found in
 * ( Chobotov, 2000 ).
 * \param meanAnomaly Mean anomaly.
 * \param pointerToCentralBody Pointer to central body.
 * \param semiMajorAxis Semi-major axis.
 * \return Elapsed time.
 */
double convertMeanAnomalyToElapsedTimeForEllipticalOrbits(
        const double& meanAnomaly, CelestialBody* pointerToCentralBody,
        const double& semiMajorAxis );

//! Convert elapsed time to mean anomaly for hyperbolic orbits.
/*!
 * Converts elapsed time to mean anomaly for hyperbolic orbits
 * ( eccentricity > 1.0 ). The equation used can be found in
 * ( Chobotov, 2000 ).
 * \param elapsedTime Elapsed time.
 * \param pointerToCentralBody Pointer to central body.
 * \param semiMajorAxis Semi-major axis.
 * \return Mean anomaly.
 */
double convertElapsedTimeToMeanAnomalyForHyperbolicOrbits(
        const double& elapsedTime, CelestialBody* pointerToCentralBody,
        const double& semiMajorAxis );

//! Convert mean anomaly to elapsed time for hyperbolic orbits.
/*!
 * Converts mean anomaly to elapsed time for hyperbolic orbits
 * ( eccentricity > 1.0 ). The equation used can be found in
 * ( Chobotov, 2000 ).
 * \param meanAnomaly Mean anomaly.
 * \param pointerToCentralBody Pointer to central body.
 * \param semiMajorAxis Semi-major axis.
 * \return Elapsed time.
 */
double convertMeanAnomalyToElapsedTimeForHyperbolicOrbits(
        const double& meanAnomaly, CelestialBody* pointerToCentralBody,
        const double& semiMajorAxis );

//! Convert mean motion to semi-major axis.
/*!
 * Converts mean motion to semi-major axis.
 * \param meanMotion Mean motion.
 * \param pointerToCentralBody Pointer to central body.
 * \return semiMajorAxis Semi-major axis.
 */
double convertMeanMotionToSemiMajorAxis( const double& meanMotion,
                                         CelestialBody* pointerToCentralBody );
}

#endif // ORBITALELEMENTCONVERSIONS_H

// End of file.

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

#include <boost/function.hpp>
#include <boost/math/special_functions/atanh.hpp>

#include <cmath>
#include <limits>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{

namespace orbital_element_conversions
{

//! Function to compute orbit semi-latus rectum.
/*!
 * Function to compute orbit semi-latus rectum.
 * \param eccentricity Eccentricity of orbit
 * \param semiMajorAxis Semi-major axis of orbit (in Tudat, this input must equal semi-latus rectum for parabolic orbits)
 * \param tolerance Eccentricity tolerance for which orbit is deemed to be parabolic.
 * \return Orbit semi-latus rectum
 */
template< typename ScalarType = double >
ScalarType computeSemiLatusRectum(
        const ScalarType eccentricity,
        const ScalarType semiMajorAxis,
        const ScalarType tolerance )
{
    // Declare semi-latus rectum.
    ScalarType semiLatusRectum;

    // Compute semi-latus rectum in the case it is not a parabola.
    if ( std::fabs( eccentricity - mathematical_constants::getFloatingInteger< ScalarType >( 1 ) ) >
         tolerance  )
    {
        semiLatusRectum = semiMajorAxis * (
                    mathematical_constants::getFloatingInteger< ScalarType >( 1 ) -
                    eccentricity * eccentricity );
    }

    // Else set the semi-latus rectum given for a parabola as the first element in the vector
    // of Keplerian elements..
    else
    {
        semiLatusRectum = semiMajorAxis;
    }
    return semiLatusRectum;
}

//! Function to compute orbit angular momentum per unit mass
/*!
 * Function to compute orbit angular momentum per unit mass
 * \param semiLatusRectum Semi-latus rectum of orbit
 * \param centralBodyGravitationalParameter Gravitational parameter of central body
 * \return Orbital angular momentum
 */
template< typename ScalarType = double >
ScalarType computeOrbitalAngularMomentumPerUnitMass(
        const ScalarType semiLatusRectum,
        const ScalarType centralBodyGravitationalParameter )
{
    return std::sqrt( semiLatusRectum * centralBodyGravitationalParameter );
}

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
template< typename ScalarType = double >
Eigen::Matrix< ScalarType, 6, 1 > convertKeplerianToCartesianElements(
        const Eigen::Matrix< ScalarType, 6, 1 >& keplerianElements,
        const ScalarType centralBodyGravitationalParameter )
{
    using std::cos;
    using std::fabs;
    using std::pow;
    using std::sin;
    using std::sqrt;

    // Set tolerance.
    const ScalarType tolerance_ = std::numeric_limits< ScalarType >::epsilon( );

    // Set local keplerian elements.
    ScalarType semiMajorAxis_ = keplerianElements( semiMajorAxisIndex );
    ScalarType eccentricity_ = keplerianElements( eccentricityIndex );
    ScalarType inclination_ = keplerianElements( inclinationIndex );
    ScalarType argumentOfPeriapsis_ = keplerianElements( argumentOfPeriapsisIndex );
    ScalarType longitudeOfAscendingNode_ = keplerianElements( longitudeOfAscendingNodeIndex );
    ScalarType trueAnomaly_ = keplerianElements( trueAnomalyIndex );

    // Pre-compute sines and cosines of involved angles for efficient computation.
    ScalarType cosineOfInclination_ = cos( inclination_ );
    ScalarType sineOfInclination_ = sin( inclination_ );
    ScalarType cosineOfArgumentOfPeriapsis_ = cos( argumentOfPeriapsis_ );
    ScalarType sineOfArgumentOfPeriapsis_ = sin( argumentOfPeriapsis_ );
    ScalarType cosineOfLongitudeOfAscendingNode_ = cos( longitudeOfAscendingNode_ );
    ScalarType sineOfLongitudeOfAscendingNode_ = sin( longitudeOfAscendingNode_ );
    ScalarType cosineOfTrueAnomaly_ = cos( trueAnomaly_ );
    ScalarType sineOfTrueAnomaly_ = sin( trueAnomaly_ );

    // Declare semi-latus rectum.
    ScalarType semiLatusRectum_ = computeSemiLatusRectum< ScalarType >( eccentricity_, semiMajorAxis_, tolerance_ );

    // Definition of position in the perifocal coordinate system.
    Eigen::Matrix< ScalarType, 2, 1 > positionPerifocal_;
    positionPerifocal_.x( ) = semiLatusRectum_ * cosineOfTrueAnomaly_
            / ( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) +
                eccentricity_ * cosineOfTrueAnomaly_ );
    positionPerifocal_.y( ) = semiLatusRectum_ * sineOfTrueAnomaly_
            / ( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) +
                eccentricity_ * cosineOfTrueAnomaly_ );

    // Definition of velocity in the perifocal coordinate system.
    Eigen::Matrix< ScalarType, 2, 1 > velocityPerifocal_(
                -sqrt( centralBodyGravitationalParameter / semiLatusRectum_ ) * sineOfTrueAnomaly_,
                sqrt( centralBodyGravitationalParameter / semiLatusRectum_ )
                * ( eccentricity_ + cosineOfTrueAnomaly_ ) );

    // Definition of the transformation matrix.
    Eigen::Matrix< ScalarType, 3, 2 > transformationMatrix_;

    // Compute the transformation matrix.
    transformationMatrix_( 0, 0 ) = cosineOfLongitudeOfAscendingNode_
            * cosineOfArgumentOfPeriapsis_ -sineOfLongitudeOfAscendingNode_
            * sineOfArgumentOfPeriapsis_ * cosineOfInclination_;
    transformationMatrix_( 0, 1 ) = -cosineOfLongitudeOfAscendingNode_
            * sineOfArgumentOfPeriapsis_ -sineOfLongitudeOfAscendingNode_
            * cosineOfArgumentOfPeriapsis_ * cosineOfInclination_;
    transformationMatrix_( 1, 0 ) = sineOfLongitudeOfAscendingNode_
            * cosineOfArgumentOfPeriapsis_ + cosineOfLongitudeOfAscendingNode_
            * sineOfArgumentOfPeriapsis_ * cosineOfInclination_;
    transformationMatrix_( 1, 1 ) = -sineOfLongitudeOfAscendingNode_
            * sineOfArgumentOfPeriapsis_ + cosineOfLongitudeOfAscendingNode_
            * cosineOfArgumentOfPeriapsis_ * cosineOfInclination_;
    transformationMatrix_( 2, 0 ) = sineOfArgumentOfPeriapsis_ * sineOfInclination_;
    transformationMatrix_( 2, 1 ) = cosineOfArgumentOfPeriapsis_ * sineOfInclination_;

    // Declare converted Cartesian elements.
    Eigen::Matrix< ScalarType, 6, 1 > convertedCartesianElements_;

    // Compute value of position in Cartesian coordinates.
    Eigen::Matrix< ScalarType, 3, 1 > position_ = transformationMatrix_ * positionPerifocal_;
    convertedCartesianElements_.segment( 0, 3 ) = position_;

    // Compute value of velocity in Cartesian coordinates.
    Eigen::Matrix< ScalarType, 3, 1 > velocity_ = transformationMatrix_ * velocityPerifocal_;
    convertedCartesianElements_.segment( 3, 3 ) = velocity_;

    // Return Cartesian elements.
    return convertedCartesianElements_;
}

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
template< typename ScalarType = double >
Eigen::Matrix< ScalarType, 6, 1 > convertCartesianToKeplerianElements(
        const Eigen::Matrix< ScalarType, 6, 1 >& cartesianElements,
        const ScalarType centralBodyGravitationalParameter )
{
    // Set tolerance.
    const ScalarType tolerance = 20.0 * std::numeric_limits< ScalarType >::epsilon( );

    // Declare converted Keplerian elements.
    Eigen::Matrix< ScalarType, 6, 1 > computedKeplerianElements_;

    // Set position and velocity vectors.
    const Eigen::Matrix< ScalarType, 3, 1 > position_( cartesianElements.segment( 0, 3 ) );
    const Eigen::Matrix< ScalarType, 3, 1 > velocity_( cartesianElements.segment( 3, 3 ) );

    // Compute orbital angular momentum vector.
    const Eigen::Matrix< ScalarType, 3, 1 > angularMomentum_( position_.cross( velocity_ ) );

    // Compute semi-latus rectum.
    const ScalarType semiLatusRectum_ = angularMomentum_.squaredNorm( )
            / centralBodyGravitationalParameter;

    // Compute unit vector to ascending node.
    Eigen::Matrix< ScalarType, 3, 1 > unitAscendingNodeVector_(
                ( Eigen::Matrix< ScalarType, 3, 1 >::UnitZ( ).cross(
                      angularMomentum_.normalized( ) ) ).normalized( ) );

    // Compute eccentricity vector.
    Eigen::Matrix< ScalarType, 3, 1 > eccentricityVector_(
                velocity_.cross( angularMomentum_ ) / centralBodyGravitationalParameter
                - position_.normalized( ) );

    // Store eccentricity.
    computedKeplerianElements_( eccentricityIndex ) = eccentricityVector_.norm( );

    // Compute and store semi-major axis.
    // Check if orbit is parabolic. If it is, store the semi-latus rectum instead of the
    // semi-major axis.
    if ( std::fabs( computedKeplerianElements_( eccentricityIndex ) -
                    mathematical_constants::getFloatingInteger< ScalarType >( 1 ) ) < tolerance )
    {
        computedKeplerianElements_( semiLatusRectumIndex ) = semiLatusRectum_;
    }

    // Else the orbit is either elliptical or hyperbolic, so store the semi-major axis.
    else
    {
        computedKeplerianElements_( semiMajorAxisIndex ) = semiLatusRectum_
                / ( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) -
                    computedKeplerianElements_( eccentricityIndex )
                    * computedKeplerianElements_( eccentricityIndex ) );
    }

    // Compute and store inclination.
    computedKeplerianElements_( inclinationIndex ) = std::acos( angularMomentum_.z( )
                                                                / angularMomentum_.norm( ) );

    // Compute and store longitude of ascending node.
    // Define the quadrant condition for the argument of perigee.
    ScalarType argumentOfPeriapsisQuandrantCondition = eccentricityVector_.z( );

    // Check if the orbit is equatorial. If it is, set the vector to the line of nodes to the
    // x-axis.
    if ( std::fabs( computedKeplerianElements_( inclinationIndex ) ) < tolerance )
    {
        unitAscendingNodeVector_ = Eigen::Matrix< ScalarType, 3, 1 >::UnitX( );

        // If the orbit is equatorial, eccentricityVector_.z( ) is zero, therefore the quadrant
        // condition is taken to be the y-component, eccentricityVector_.y( ).
        argumentOfPeriapsisQuandrantCondition = eccentricityVector_.y( );
    }

    // Compute and store the resulting longitude of ascending node.
    computedKeplerianElements_( longitudeOfAscendingNodeIndex )
            = std::acos( unitAscendingNodeVector_.x( ) );

    // Check if the quandrant is correct.
    if ( unitAscendingNodeVector_.y( ) <
         mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        computedKeplerianElements_( longitudeOfAscendingNodeIndex ) =
                mathematical_constants::getFloatingInteger< ScalarType >( 2 ) *
                mathematical_constants::getPi< ScalarType >( ) -
                computedKeplerianElements_( longitudeOfAscendingNodeIndex );
    }

    // Compute and store argument of periapsis.
    // Define the quadrant condition for the true anomaly.
    ScalarType trueAnomalyQuandrantCondition = position_.dot( velocity_ );

    // Check if the orbit is circular. If it is, set the eccentricity vector to unit vector
    // pointing to the ascending node, i.e. set the argument of periapsis to zero.
    if ( std::fabs( computedKeplerianElements_( eccentricityIndex ) ) < tolerance )
    {
        eccentricityVector_ = unitAscendingNodeVector_;

        computedKeplerianElements_( argumentOfPeriapsisIndex ) =
                mathematical_constants::getFloatingInteger< ScalarType >( 0 );

        // Check if orbit is also equatorial and set true anomaly quandrant check condition
        // accordingly.
        if ( unitAscendingNodeVector_ == Eigen::Matrix< ScalarType, 3, 1 >::UnitX( ) )
        {
            // If the orbit is circular, position_.dot( velocity_ ) = 0, therefore this value
            // cannot be used as a quadrant condition. Moreover, if the orbit is equatorial,
            // position_.z( ) is also zero and therefore the quadrant condition is taken to be the
            // y-component, position_.y( ).
            trueAnomalyQuandrantCondition = position_.y( );
        }

        else
        {
            // If the orbit is circular, position_.dot( velocity_ ) = 0, therefore the quadrant
            // condition is taken to be the z-component of the position, position_.z( ).
            trueAnomalyQuandrantCondition = position_.z( );
        }
    }

    // Else, compute the argument of periapsis as the angle between the eccentricity vector and
    // the unit vector to the ascending node.
    else
    {
        ScalarType eccentricityAscendingNodeDotProduct = eccentricityVector_.normalized( ).dot( unitAscendingNodeVector_ );

        // Check whether dot product is in bounds (might be out of bounds due to numerical noise).
        if( eccentricityAscendingNodeDotProduct < mathematical_constants::getFloatingInteger< ScalarType >( -1 ) )
        {
            computedKeplerianElements_( argumentOfPeriapsisIndex ) = mathematical_constants::getPi< ScalarType >( );
        }
        else if( eccentricityAscendingNodeDotProduct > mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
        {
            computedKeplerianElements_( argumentOfPeriapsisIndex ) =
                    mathematical_constants::getFloatingInteger< ScalarType >( 0 );
        }
        else
        {
             computedKeplerianElements_( argumentOfPeriapsisIndex ) = std::acos( eccentricityAscendingNodeDotProduct );
        }

        // Check if the quadrant is correct.
        if ( argumentOfPeriapsisQuandrantCondition <
             mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
        {
            computedKeplerianElements_( argumentOfPeriapsisIndex ) =
                    mathematical_constants::getFloatingInteger< ScalarType >( 2 ) *
                    mathematical_constants::getPi< ScalarType >( ) -
                    computedKeplerianElements_( argumentOfPeriapsisIndex );
        }
    }

    // Compute dot-product of position and eccentricity vectors.
    ScalarType dotProductPositionAndEccentricityVectors
            = position_.normalized( ).dot( eccentricityVector_.normalized( ) );

    // Check if the dot-product is one of the limiting cases: 0.0, -1.0 or 1.0
    // (within prescribed tolerance).
    if ( std::fabs( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) -
                    dotProductPositionAndEccentricityVectors ) < tolerance )
    {
        dotProductPositionAndEccentricityVectors =
                mathematical_constants::getFloatingInteger< ScalarType >( 1 );
    }

    if ( std::fabs( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) +
                    dotProductPositionAndEccentricityVectors ) < tolerance )
    {
        dotProductPositionAndEccentricityVectors =
                -mathematical_constants::getFloatingInteger< ScalarType >( 1 );
    }

    if ( std::fabs( dotProductPositionAndEccentricityVectors ) < tolerance )
    {
        dotProductPositionAndEccentricityVectors  =
                mathematical_constants::getFloatingInteger< ScalarType >( 0 );
    }

    // Compute and store true anomaly.
    computedKeplerianElements_( trueAnomalyIndex )
            = std::acos( dotProductPositionAndEccentricityVectors );

    // Check if the quandrant is correct.
    if ( trueAnomalyQuandrantCondition <
         mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        computedKeplerianElements_( trueAnomalyIndex ) =
                mathematical_constants::getFloatingInteger< ScalarType >( 2 ) *
                mathematical_constants::getPi< ScalarType >( ) -
                computedKeplerianElements_( trueAnomalyIndex );
    }

    // Return converted Keplerian elements.
    return computedKeplerianElements_;
}

//! Convert Cartesian to Keplerian orbital elements.
/*!
 * Converts Cartesian to Keplerian orbital elements, using function pointers to retrieve the cartesian state and gravitational
 * parameter.
 * \param cartesianElementsFunction Function that returns vector containing Cartesian elements.
 * \param centralBodyGravitationalParameterFunction Function that returns  gravitational parameter of central body.
 * \return Converted state in Keplerian elements.
 * \sa convertCartesianToKeplerianElements
 */
template< typename ScalarType = double >
Eigen::Matrix< ScalarType, 6, 1 > convertCartesianToKeplerianElementsFromFunctions(
        const boost::function< Eigen::Matrix< ScalarType, 6, 1 >( ) > cartesianElementsFunction,
        const boost::function< ScalarType( ) > centralBodyGravitationalParameterFunction )
{
    return convertCartesianToKeplerianElements(
                cartesianElementsFunction( ), centralBodyGravitationalParameterFunction( ) );
}

//! Convert true anomaly to (elliptical) eccentric anomaly.
/*!
 * Converts true anomaly to eccentric anomaly for elliptical orbits ( 0 <= eccentricity < 1.0 ).
 * The equations used can be found in (Chobotov, 2002).
 * \param trueAnomaly True anomaly.                                                           [rad]
 * \param eccentricity Eccentricity.                                                            [-]
 * \return (Elliptical) Eccentric anomaly.                                                    [rad]
 */
template< typename ScalarType = double >
ScalarType convertTrueAnomalyToEllipticalEccentricAnomaly(
        const ScalarType trueAnomaly, const ScalarType eccentricity )

{
    using std::cos;
    using std::sqrt;
    using std::atan2;

    if ( eccentricity >= mathematical_constants::getFloatingInteger< ScalarType >( 1 ) ||
         eccentricity < mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        throw std::runtime_error( "Eccentricity is invalid when converting true to elliptical eccentric anomaly." );
    }
    else
    {
        // Declare and compute sine and cosine of eccentric anomaly.
        ScalarType sineOfEccentricAnomaly_ =
                sqrt( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) -
                      eccentricity * eccentricity ) * std::sin( trueAnomaly ) /
                ( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) +
                  eccentricity * cos( trueAnomaly ) );
        ScalarType cosineOfEccentricAnomaly_ = ( eccentricity + cos( trueAnomaly ) )
                / ( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) +
                    eccentricity * cos( trueAnomaly ) );

        // Return elliptic eccentric anomaly.
        return atan2( sineOfEccentricAnomaly_, cosineOfEccentricAnomaly_ );
    }
}

//! Convert true anomaly to hyperbolic eccentric anomaly.
/*!
 * Converts true anomaly to hyperbolic eccentric anomaly for hyperbolic orbits
 * ( eccentricity > 1.0 ). The equations used can be found in (Chobotov, 2002).
 * \param trueAnomaly True anomaly.                                                           [rad]
 * \param eccentricity Eccentricity.                                                            [-]
 * \return Hyperbolic eccentric anomaly.                                                      [rad]
 */
template< typename ScalarType = double >
ScalarType convertTrueAnomalyToHyperbolicEccentricAnomaly( const ScalarType trueAnomaly,
                                                           const ScalarType eccentricity )
{
    if ( eccentricity <= mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
    {
        throw std::runtime_error( "Eccentricity is invalid  when converting true to hyperbolic eccentric anomaly." );
    }

    else
    {
        using std::cos;

        // Compute hyperbolic sine and hyperbolic cosine of hyperbolic eccentric anomaly.
        ScalarType hyperbolicSineOfHyperbolicEccentricAnomaly_
                = std::sqrt( eccentricity * eccentricity -
                             mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
                * std::sin( trueAnomaly ) /
                ( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) +
                  cos( trueAnomaly ) );

        ScalarType hyperbolicCosineOfHyperbolicEccentricAnomaly_
                = ( cos( trueAnomaly ) + eccentricity ) /
                ( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) +
                  cos( trueAnomaly ) );

        // Return hyperbolic eccentric anomaly.
        return boost::math::atanh( hyperbolicSineOfHyperbolicEccentricAnomaly_
                                   / hyperbolicCosineOfHyperbolicEccentricAnomaly_ );
    }
}

//! Convert true anomaly to eccentric anomaly.
/*!
 * Converts true anomaly to eccentric anomaly for elliptical and hyperbolic orbits
 * ( eccentricity < 1.0 && eccentricity > 1.0 ). This function is essentially a wrapper for
 * convertTrueAnomalyToEllipticalEccentricAnomaly() and
 * convertTrueAnomalyToHyperbolicEccentricAnomaly(). It should be used in cases where the
 * eccentricity of the orbit is not known a priori. Currently, this implementation performs a
 * check on the eccentricity and throws an error for eccentricity < 0.0 and parabolic orbits, which
 * have not been implemented. The equations used can be found in (Chobotov, 2002).
 * \param trueAnomaly True anomaly.                                                           [rad]
 * \param eccentricity Eccentricity.                                                            [-]
 * \return Eccentric anomaly.                                                                 [rad]
 */
template< typename ScalarType = double >
ScalarType convertTrueAnomalyToEccentricAnomaly( const ScalarType trueAnomaly,
                                                 const ScalarType eccentricity )
{
    // Declare computed eccentric anomaly.
    ScalarType eccentricAnomaly_ = 0.0;

    // Check if eccentricity is invalid and throw an error if true.
    if ( eccentricity < mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        throw std::runtime_error( "Eccentricity is invalid when converting true to eccentric anomaly." );
    }

    // Check if orbit is parabolic and throw an error if true.
    else if ( std::fabs( eccentricity -
                         mathematical_constants::getFloatingInteger< ScalarType >( 1 ) ) <
              std::numeric_limits< ScalarType >::epsilon( ) )
    {
        throw std::runtime_error( "Parabolic orbits have not yet been implemented when converting true to eccentric anomaly." );
    }

    // Check if orbit is elliptical and compute eccentric anomaly.
    else if ( eccentricity >= mathematical_constants::getFloatingInteger< ScalarType >( 0 ) &&
              eccentricity < mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
    {
        eccentricAnomaly_ = convertTrueAnomalyToEllipticalEccentricAnomaly< ScalarType >(
                    trueAnomaly, eccentricity );
    }

    else if ( eccentricity > mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
    {
        eccentricAnomaly_ = convertTrueAnomalyToHyperbolicEccentricAnomaly< ScalarType >(
                    trueAnomaly, eccentricity );
    }

    // Return computed eccentric anomaly.
    return eccentricAnomaly_;
}

//! Convert (elliptical) eccentric anomaly to true anomaly.
/*!
 * Converts eccentric anomaly to true anomaly for elliptical orbits ( 0 <= eccentricity < 1.0 ).
 * The equations used can be found in (Chobotov, 2002).
 * \param ellipticEccentricAnomaly Elliptical eccentric anomaly.                              [rad]
 * \param eccentricity Eccentricity.                                                            [-]
 * \return True anomaly.                                                                      [rad]
 */
template< typename ScalarType = double >
ScalarType convertEllipticalEccentricAnomalyToTrueAnomaly(
        const ScalarType ellipticEccentricAnomaly,
        const ScalarType eccentricity )
{
    if ( eccentricity >= mathematical_constants::getFloatingInteger< ScalarType >( 1 ) ||
         eccentricity < mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        throw std::runtime_error( "Eccentricity is invalid when converting elliptical eccentric to true anomaly." );
    }

    else
    {
        using std::cos;
        using std::sqrt;

        // Compute sine and cosine of true anomaly.
        ScalarType sineOfTrueAnomaly_ =
                sqrt( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) -
                      eccentricity * eccentricity ) *
                std::sin( ellipticEccentricAnomaly )
                / ( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) -
                    eccentricity * cos( ellipticEccentricAnomaly ) );

        ScalarType cosineOfTrueAnomaly_ = ( cos( ellipticEccentricAnomaly ) - eccentricity )
                / ( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) -
                    eccentricity * cos( ellipticEccentricAnomaly ) );

        // Return true anomaly.
        return std::atan2( sineOfTrueAnomaly_, cosineOfTrueAnomaly_  );
    }
}

//! Convert hyperbolic eccentric anomaly to true anomaly.
/*!
 * Converts hyperbolic eccentric anomaly to true anomaly for hyperbolic orbits
 * ( eccentricity > 1.0 ). The equations used can be found in (Chobotov, 2002).
 * \param hyperbolicEccentricAnomaly Hyperbolic eccentric anomaly.                            [rad]
 * \param eccentricity Eccentricity.                                                            [-]
 * \return True anomaly.                                                                      [rad]
 */
template< typename ScalarType = double >
ScalarType convertHyperbolicEccentricAnomalyToTrueAnomaly(
        const ScalarType hyperbolicEccentricAnomaly,
        const ScalarType eccentricity )
{
    if ( eccentricity <= mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
    {
         throw std::runtime_error( "Eccentricity is invalid  when converting hyperbolic eccentric to true anomaly." );
    }

    else
    {
        using std::cosh;

        // Compute sine and cosine of true anomaly.
        ScalarType sineOfTrueAnomaly_
                = std::sqrt( eccentricity * eccentricity -
                             mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
                * std::sinh( hyperbolicEccentricAnomaly )
                / ( eccentricity * cosh( hyperbolicEccentricAnomaly ) -
                    mathematical_constants::getFloatingInteger< ScalarType >( 1 ) );

        ScalarType cosineOfTrueAnomaly_
                = ( eccentricity - cosh( hyperbolicEccentricAnomaly ) )
                / ( eccentricity * cosh( hyperbolicEccentricAnomaly ) -
                    mathematical_constants::getFloatingInteger< ScalarType >( 1 ) );

        // Return true anomaly.
        return std::atan2( sineOfTrueAnomaly_, cosineOfTrueAnomaly_ );
    }

}


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
template< typename ScalarType = double >
ScalarType convertEccentricAnomalyToTrueAnomaly( const ScalarType eccentricAnomaly,
                                                 const ScalarType eccentricity )
{
    // Declare computed true anomaly.
    ScalarType trueAnomaly_ = -mathematical_constants::getFloatingInteger< ScalarType >( 0 );

    // Check if eccentricity is invalid and throw an error if true.
    if ( eccentricity < mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        throw std::runtime_error( "Eccentricity is invalid when converting eccentric to true anomaly." );
    }

    // Check if orbit is parabolic and throw an error if true.
    else if ( std::fabs( eccentricity -
                         mathematical_constants::getFloatingInteger< ScalarType >( 1 ) ) <
              std::numeric_limits< ScalarType >::epsilon( ) )
    {
        throw std::runtime_error( "Parabolic orbits have not yet been implemented when converting eccentric to true anomaly." );
    }

    // Check if orbit is elliptical and compute true anomaly.
    else if ( eccentricity >= mathematical_constants::getFloatingInteger< ScalarType >( 0 ) &&
              eccentricity < mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
    {
        trueAnomaly_ = convertEllipticalEccentricAnomalyToTrueAnomaly( eccentricAnomaly,
                                                                       eccentricity );
    }

    else if ( eccentricity > mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
    {
        trueAnomaly_ = convertHyperbolicEccentricAnomalyToTrueAnomaly( eccentricAnomaly,
                                                                       eccentricity );
    }

    // Return computed true anomaly.
    return trueAnomaly_;
}


//! Convert (elliptical) eccentric anomaly to mean anomaly.
/*!
 * Converts eccentric anomaly to mean anomaly for elliptical orbits ( 0 <= eccentricity < 1.0 ).
 * The equations used can be found in (Chobotov, 2002).
 * \param eccentricity Eccentricity.                                                            [-]
 * \param ellipticalEccentricAnomaly (Elliptical) eccentric anomaly [rad].
 * \return Mean anomaly [rad].
 */
template< typename ScalarType = double >
ScalarType convertEllipticalEccentricAnomalyToMeanAnomaly(
        const ScalarType ellipticalEccentricAnomaly,
        const ScalarType eccentricity )
{
    return ellipticalEccentricAnomaly - eccentricity * std::sin( ellipticalEccentricAnomaly );
}


//! Convert hyperbolic eccentric anomaly to mean anomaly.
/*!
 * Converts hyperbolic eccentric anomaly to mean anomaly for hyperbolic orbits
 * ( eccentricity > 1.0 ). The equations used can be found in (Chobotov, 2002).
 * \param hyperbolicEccentricAnomaly Hyperbolic eccentric anomaly.                            [rad]
 * \param eccentricity Eccentricity.                                                            [-]
 * \return Mean anomaly.                                                                      [rad]
 */
template< typename ScalarType = double >
ScalarType convertHyperbolicEccentricAnomalyToMeanAnomaly(
        const ScalarType hyperbolicEccentricAnomaly,
        const ScalarType eccentricity )
{
    return eccentricity * std::sinh( hyperbolicEccentricAnomaly ) - hyperbolicEccentricAnomaly;
}

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
template< typename ScalarType = double >
ScalarType convertEccentricAnomalyToMeanAnomaly(
        const ScalarType eccentricAnomaly,
        const ScalarType eccentricity )
{
    // Declare computed mean anomaly.
    ScalarType meanAnomaly_ = 0.0;

    // Check if eccentricity is invalid and throw an error if true.
    if ( eccentricity < mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        throw std::runtime_error( "Eccentricity is invalid when converting eccentric to mean anomaly." );
    }

    // Check if orbit is parabolic and throw an error if true.
    else if ( std::fabs( eccentricity -
                         mathematical_constants::getFloatingInteger< ScalarType >( 1 ) ) <
              std::numeric_limits< ScalarType >::epsilon( ) )
    {
        throw std::runtime_error(
                            "Parabolic orbits have not yet been implemented when converting eccentric to mean anomaly." );
    }

    // Check if orbit is elliptical and compute true anomaly.
    else if ( eccentricity >=
              mathematical_constants::getFloatingInteger< ScalarType >( 0 ) &&
              eccentricity < mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
    {
        meanAnomaly_ = convertEllipticalEccentricAnomalyToMeanAnomaly< ScalarType >(
                    eccentricAnomaly, eccentricity );
    }

    else if ( eccentricity > mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
    {
        meanAnomaly_ = convertHyperbolicEccentricAnomalyToMeanAnomaly< ScalarType >(
                    eccentricAnomaly, eccentricity );
    }

    // Return computed mean anomaly.
    return meanAnomaly_;
}

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
template< typename ScalarType = double >
ScalarType convertElapsedTimeToEllipticalMeanAnomalyChange(
        const ScalarType elapsedTime, const ScalarType centralBodyGravitationalParameter,
        const ScalarType semiMajorAxis )
{
    // Check if semi-major axis is invalid and throw error if true.
    if ( semiMajorAxis < mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
       throw std::runtime_error( "Semi-major axis is invalid when converting elapsed time to mean anomlay change to." );
    }

    // Else return elliptical mean anomaly change.
    else
    {
        return std::sqrt( centralBodyGravitationalParameter
                          / ( semiMajorAxis * semiMajorAxis * semiMajorAxis ) ) * elapsedTime;
    }
}


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
template< typename ScalarType = double >
ScalarType convertElapsedTimeToHyperbolicMeanAnomalyChange(
        const ScalarType elapsedTime, const ScalarType centralBodyGravitationalParameter,
        const ScalarType semiMajorAxis )
{
    // Check if semi-major axis is invalid and throw error if true.
    if ( semiMajorAxis > mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        throw std::runtime_error(
                    "Semi-major axis is invalid when converting elapsed time to hyperbolic mean anomlay change to." );
    }

    // Else return hyperbolic mean anomaly change.
    else
    {
        return std::sqrt( centralBodyGravitationalParameter
                          / ( - semiMajorAxis * semiMajorAxis * semiMajorAxis ) ) * elapsedTime;
    }
}

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
template< typename ScalarType = double >
ScalarType convertElapsedTimeToMeanAnomalyChange(
        const ScalarType elapsedTime, const ScalarType centralBodyGravitationalParameter,
        const ScalarType semiMajorAxis )
{
    // Declare computed mean anomaly change.
    ScalarType meanAnomalyChange_ = -mathematical_constants::getFloatingInteger< ScalarType >( 0 );

    // Check if orbit is elliptical and compute mean anomaly change.
    if ( semiMajorAxis > mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        meanAnomalyChange_ = convertElapsedTimeToEllipticalMeanAnomalyChange(
                    elapsedTime, centralBodyGravitationalParameter, semiMajorAxis );
    }

    // Else orbit is hyperbolic; compute mean anomaly change.
    else if ( semiMajorAxis < mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        meanAnomalyChange_ = convertElapsedTimeToHyperbolicMeanAnomalyChange(
                    elapsedTime, centralBodyGravitationalParameter, semiMajorAxis );
    }

    // Return computed mean anomaly change.
    return meanAnomalyChange_;
}


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
template< typename ScalarType = double >
ScalarType convertEllipticalMeanAnomalyChangeToElapsedTime(
        const ScalarType ellipticalMeanAnomalyChange,
        const ScalarType centralBodyGravitationalParameter,
        const ScalarType semiMajorAxis )
{
    // Check if semi-major axis is invalid and throw error if true.
    if ( semiMajorAxis < mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        throw std::runtime_error( "Semi-major axis is invalid when converting mean anomlay change to elapsed time."  );
    }

    // Else return elapsed time.
    else
    {
        return ellipticalMeanAnomalyChange * std::sqrt(
                    semiMajorAxis * semiMajorAxis * semiMajorAxis
                    / centralBodyGravitationalParameter );
    }
}

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
template< typename ScalarType = double >
ScalarType convertHyperbolicMeanAnomalyChangeToElapsedTime(
        const ScalarType hyperbolicMeanAnomalyChange,
        const ScalarType centralBodyGravitationalParameter,
        const ScalarType semiMajorAxis )
{
    // Check if semi-major axis is invalid and throw error if true.
    if ( semiMajorAxis > mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
       throw std::runtime_error( "Semi-major axis is invalid when converting hyperbolic mean anomlay change to elapsed time." );
    }

    // Else return elapsed time.
    else
    {
        return std::sqrt( -semiMajorAxis * semiMajorAxis * semiMajorAxis
                          / centralBodyGravitationalParameter ) * hyperbolicMeanAnomalyChange;
    }
}

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
template< typename ScalarType = double >
ScalarType convertMeanAnomalyChangeToElapsedTime(
        const ScalarType meanAnomalyChange, const ScalarType centralBodyGravitationalParameter,
        const ScalarType semiMajorAxis )
{
    // Declare computed elapsed time.
    ScalarType elapsedTime_ = mathematical_constants::getFloatingInteger< ScalarType >( 0 );

    // Check if orbit is elliptical and compute elapsed time.
    if ( semiMajorAxis > mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        elapsedTime_ = convertEllipticalMeanAnomalyChangeToElapsedTime< ScalarType >(
                    meanAnomalyChange, centralBodyGravitationalParameter, semiMajorAxis );
    }

    // Else orbit is hyperbolic; compute elapsed time.
    else if ( semiMajorAxis < mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        elapsedTime_ = convertHyperbolicMeanAnomalyChangeToElapsedTime< ScalarType >(
                    meanAnomalyChange, centralBodyGravitationalParameter, semiMajorAxis );
    }

    // Return computed elapsed time.
    return elapsedTime_;
}

//! Convert (elliptical) mean motion to semi-major axis.
/*!
 * Converts mean motion to semi-major axis for elliptical orbits.
 * \param ellipticalMeanMotion (Elliptical) Mean motion.                                    [rad/s]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \return semiMajorAxis Semi-major axis.                                                       [m]
 */
template< typename ScalarType = double >
ScalarType convertEllipticalMeanMotionToSemiMajorAxis(
        const ScalarType ellipticalMeanMotion, const ScalarType centralBodyGravitationalParameter )
{
    return std::pow( centralBodyGravitationalParameter
                     / ( ellipticalMeanMotion * ellipticalMeanMotion ),
                     mathematical_constants::getFloatingFraction< ScalarType >( 1, 3 ) );
}

//! Convert semi-major axis to elliptical mean motion.
/*!
 * Converts semi-major axis to elliptical mean motion.
 * \param semiMajorAxis Semi-major axis.                                                        [m]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \return ellipticalMeanMotion (Elliptical) Mean motion.                                   [rad/s]
 */
template< typename ScalarType = double >
ScalarType convertSemiMajorAxisToEllipticalMeanMotion(
        const ScalarType semiMajorAxis, const ScalarType centralBodyGravitationalParameter )
{
    // Check if semi-major axis is invalid and throw error if true.
    if ( semiMajorAxis < mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        throw std::runtime_error( "Semi-major axis is invalid when converting semi major axis to elliptical mean motion." );
    }

    // Else compute and return elliptical mean motion.
    {
        return std::sqrt( centralBodyGravitationalParameter /
                          ( semiMajorAxis * semiMajorAxis * semiMajorAxis ) );
    }
}
} // namespace orbital_element_conversions

} // namespace tudat

#endif // TUDAT_ORBITAL_ELEMENT_CONVERSIONS_H

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
 *      Verified Interval Propagation, Bart Rmgens; Delft (2012?). Code archive.
 *      Modified Equinoctial Orbital Elements, author unknown;
 *          http://www.cdeagle.com/pdf/mee.pdf (2010?).
 *      Survey of Orbital Element Sets, Gerald R. Hintz; Journal of Guidance, Control and
 *          Dynamics (2008, Vol. 31 - Nr. 3).
 *      Code archive, E. Heeren (fellow Tudat developer).
 *
 */

#ifndef TUDAT_MODIFIED_EQUINOCTIAL_ELEMENT_CONVERSIONS_H
#define TUDAT_MODIFIED_EQUINOCTIAL_ELEMENT_CONVERSIONS_H

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/missionGeometry.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"

namespace tudat
{

namespace orbital_element_conversions
{


//! Convert Keplerian to modified equinoctial orbital elements using MEE explicit equation set.
/*!
 * Converts Keplerian to modified equinoctial elements using the prograde/retrograde equation
 * determined by the values of the input Kepler elements. If input exceeds allowable ranges,
 * an error is thrown. NOTE: This function automatically selects the boolean parameter that defines
 * the location of the singularity, depending on whether the orbit is prograde or retrograde.
 * \param keplerianElements Vector containing Keplerian elements. Order of elements is important!
 *         keplerianElements( 0 ) = semi-major axis,                                            [m]
 *         keplerianElements( 1 ) = eccentricity,                                               [-]
 *         keplerianElements( 2 ) = inclination (in range [0,PI]),                            [rad]
 *         keplerianElements( 3 ) = argument of periapsis,                                    [rad]
 *         keplerianElements( 4 ) = longitude of ascending node,                              [rad]
 *         keplerianElements( 5 ) = true anomaly.                                             [rad]
 * \param flipSingularityToZeroInclination Boolean flag to indicate whether the set of equations for
 *          the inclination = 180 degrees (false) or 0 degrees (true) singular case are to be used.
 *          Take note: the same set of equations
 *          is required for conversion back to Keplerian elements to retrieve original state!
 * \return Converted state in modified equinoctial elements. The order of elements is fixed!
 *         modifiedEquinoctialElements( 0 ) = semi-latus rectum,                                [m]
 *         modifiedEquinoctialElements( 1 ) = f-element,                                        [-]
 *         modifiedEquinoctialElements( 2 ) = g-element,                                        [-]
 *         modifiedEquinoctialElements( 3 ) = h-element,                                        [-]
 *         modifiedEquinoctialElements( 4 ) = k-element,                                        [-]
 *         modifiedEquinoctialElements( 5 ) = true longitude.                                 [rad]
 */
template< typename ScalarType = double >
Eigen::Matrix< ScalarType, 6, 1 > convertKeplerianToModifiedEquinoctialElements(
        const Eigen::Matrix< ScalarType, 6, 1 >& keplerianElements,
        const bool flipSingularityToZeroInclination )
{
    using mathematical_constants::getPi;
    using mathematical_constants::getFloatingInteger;

    // Declaring eventual output vector.
    Eigen::Matrix< ScalarType, 6, 1 > modifiedEquinoctialState( 6 );

    // Compute semi-latus rectum.
    ScalarType singularityTolerance = 5.0 * std::numeric_limits< ScalarType >::epsilon( );

    // Extracting eccentricity for ease of referencing.
    ScalarType eccentricity = keplerianElements( eccentricityIndex );

    // If e is (very near) one, then semi-major axis is undefined and thus the first kepler
    // element is the semi-latus rectum.
    if ( std::fabs( eccentricity - getFloatingInteger< ScalarType >( 1 ) ) < singularityTolerance )
    {
        modifiedEquinoctialState( semiLatusRectumIndex )
                = keplerianElements( semiLatusRectumIndex );
    }
    else // eccentricity is significantly away from singularity and can be computed with p=a(1-e).
    {
        modifiedEquinoctialState( semiLatusRectumIndex ) =
                keplerianElements( semiMajorAxisIndex )
                * ( getFloatingInteger< ScalarType >( 1 ) - eccentricity * eccentricity );
    }

    // Determine prograde-ness, as other five parameters depend on that fact.
    ScalarType inclination = keplerianElements( inclinationIndex );

    // If inclination is outside range [0,PI].
    if ( ( inclination < getFloatingInteger< ScalarType >( 0 ) ) || ( inclination > getPi< ScalarType >( ) ) )
    {
        // Define the error message.
        throw std::runtime_error(  "Inclination is expected in range [0," +  std::to_string( getPi< ScalarType >( ) ) + "]\n"
                     + "Specified inclination: " + std::to_string( inclination ) + " rad." );
    }
    //Else, nothing wrong and continue.

    // Compute set dependant helper parameters.
    ScalarType argumentOfPeriapsisAndAscendingNode = getFloatingInteger< ScalarType >( 0 );
    ScalarType tangentOfHalfInclination = getFloatingInteger< ScalarType >( 0 );

    // If normal set of equations is to be used (i.e. not inverse set for singularity).
    if ( !flipSingularityToZeroInclination )
    {
        // Add (+) argument of periapsis and longitude of ascending node.
        argumentOfPeriapsisAndAscendingNode = keplerianElements( argumentOfPeriapsisIndex )
                + keplerianElements( longitudeOfAscendingNodeIndex );

        // Take the tangent of half the inclination.
        tangentOfHalfInclination = std::tan( inclination / getFloatingInteger< ScalarType >( 2 ) );
    }
    else // the orbit is retrograde.
    {
        // Subtract (-) longitude of ascending node from argument of periapsis.
        argumentOfPeriapsisAndAscendingNode = keplerianElements( argumentOfPeriapsisIndex )
                - keplerianElements( longitudeOfAscendingNodeIndex );

        // Take the inverse of the tangent of half the inclination to avoid singularity at tan
        // PI/2.
        tangentOfHalfInclination = getFloatingInteger< ScalarType >( 1 ) /
                std::tan( inclination / getFloatingInteger< ScalarType >( 2 ) );
    }

    // Compute f-element.
    modifiedEquinoctialState( fElementIndex ) = eccentricity
            * std::cos( argumentOfPeriapsisAndAscendingNode );

    // Compute g-element.
    modifiedEquinoctialState( gElementIndex ) = eccentricity
            * std::sin( argumentOfPeriapsisAndAscendingNode );

    // Compute h-element.
    modifiedEquinoctialState( hElementIndex ) = tangentOfHalfInclination
            * std::cos( keplerianElements( longitudeOfAscendingNodeIndex ) );

    // Compute k-element.
    modifiedEquinoctialState( kElementIndex ) = tangentOfHalfInclination
            * std::sin( keplerianElements( longitudeOfAscendingNodeIndex ) );

    // Compute true longitude (modulo 2 PI to keep within interval -2PI to 2PI).
    modifiedEquinoctialState( trueLongitudeIndex )
            = basic_mathematics::computeModulo(
                argumentOfPeriapsisAndAscendingNode
                + keplerianElements( trueAnomalyIndex ), getFloatingInteger< ScalarType >( 2 ) * getPi< ScalarType >( ) );

    // Give back result.
    return modifiedEquinoctialState;
}

//! Convert Keplerian to modified equinoctial orbital elements using implicit MEE equation set.
/*!
 * Converts Keplerian to modified equinoctial elements using the prograde/retrograde equation
 * determined by the values of the input Kepler elements. If input exceeds allowable ranges,
 * an error is thrown.
 * \param keplerianElements Vector containing Keplerian elements. Order of elements is important!
 *         keplerianElements( 0 ) = semi-major axis,                                            [m]
 *         keplerianElements( 1 ) = eccentricity,                                               [-]
 *         keplerianElements( 2 ) = inclination (in range [0,PI]),                            [rad]
 *         keplerianElements( 3 ) = argument of periapsis,                                    [rad]
 *         keplerianElements( 4 ) = longitude of ascending node,                              [rad]
 *         keplerianElements( 5 ) = true anomaly.                                             [rad]
 * \return Converted state in modified equinoctial elements. The order of elements is fixed!
 *         modifiedEquinoctialElements( 0 ) = semi-latus rectum,                                [m]
 *         modifiedEquinoctialElements( 1 ) = f-element,                                        [-]
 *         modifiedEquinoctialElements( 2 ) = g-element,                                        [-]
 *         modifiedEquinoctialElements( 3 ) = h-element,                                        [-]
 *         modifiedEquinoctialElements( 4 ) = k-element,                                        [-]
 *         modifiedEquinoctialElements( 5 ) = true longitude.                                 [rad]
 */
template< typename ScalarType = double >
Eigen::Matrix< ScalarType, 6, 1 > convertKeplerianToModifiedEquinoctialElements(
        const Eigen::Matrix< ScalarType, 6, 1 >& keplerianElements )
// Based on Hintz, 2008.
{
    // Check if orbit is retrograde
    bool flipSingularityToZeroInclination =
            mission_geometry::isOrbitRetrograde( keplerianElements );

    // Convert to modified equinoctial elements.
    return convertKeplerianToModifiedEquinoctialElements( keplerianElements,
                                                          flipSingularityToZeroInclination );
}

//! Convert modified equinoctial to Keplerian orbital elements.
/*!
 * Converts modified equinoctial elements to Keplerian using one of two sets of equations specified
 * by the user.
 * \param modifiedEquinoctialElements Vector containing modified equinoctial elements. Order of
 *          elements is important!
 *         modifiedEquinoctialElements( 0 ) = semi-latus rectum,                                [m]
 *         modifiedEquinoctialElements( 1 ) = f-element,                                        [-]
 *         modifiedEquinoctialElements( 2 ) = g-element,                                        [-]
 *         modifiedEquinoctialElements( 3 ) = h-element,                                        [-]
 *         modifiedEquinoctialElements( 4 ) = k-element,                                        [-]
 *         modifiedEquinoctialElements( 5 ) = true longitude.                                 [rad]
 * \param flipSingularityToZeroInclination Boolean flag to indicate whether the set of equations for
 *          the inclination = 180 degrees (false) or 0 degrees (true) singular case are to be used.
 *          Take note: the same set of equations
 *          is required for conversion back to modified equinoctial elements to retrieve original
 *          state!
 * \return Converted state in Keplerian elements. The order of elements is fixed!
 *         keplerianElements( 0 ) = semi-major axis,                                            [m]
 *         keplerianElements( 1 ) = eccentricity,                                               [-]
 *         keplerianElements( 2 ) = inclination [0, 180],                                     [rad]
 *         keplerianElements( 3 ) = argument of periapsis,                                    [rad]
 *         keplerianElements( 4 ) = longitude of ascending node,                              [rad]
 *         keplerianElements( 5 ) = true anomaly.                                             [rad]
 * \note Using unknown source pdf, code archive E. Heeren and personal derivation based on Hintz 2008.
 */
template< typename ScalarType = double >
Eigen::Matrix< ScalarType, 6, 1 > convertModifiedEquinoctialToKeplerianElements(
        const Eigen::Matrix< ScalarType, 6, 1 >& modifiedEquinoctialElements,
        const bool flipSingularityToZeroInclination )
{
    using mathematical_constants::getPi;
    using mathematical_constants::getFloatingInteger;

    // Declaration of output vector.
    Eigen::Matrix< ScalarType, 6, 1 > convertedKeplerianElements = Eigen::Matrix< ScalarType, 6, 1 >::Zero( 6 );

    // for ease of referencing, almost all modified equinoctial elements.
    ScalarType fElement = modifiedEquinoctialElements( fElementIndex );
    ScalarType gElement = modifiedEquinoctialElements( gElementIndex );
    ScalarType hElement = modifiedEquinoctialElements( hElementIndex );
    ScalarType kElement = modifiedEquinoctialElements( kElementIndex );

    // Tolerance for singularities.
    ScalarType singularityTolerance = 5.0 * std::numeric_limits< ScalarType >::epsilon( );

    // Compute eccentricity.0
    ScalarType eccentricity = std::sqrt( fElement * fElement + gElement * gElement );
    convertedKeplerianElements( eccentricityIndex ) = eccentricity;

    // Compute semi-major axis.
    // If eccentricity is not near-parabolic.
    if ( std::fabs( eccentricity - getFloatingInteger< ScalarType >( 1 ) ) > singularityTolerance )
    {
        // Use semi-latus rectum and eccentricity to calculate semi-major axis with a=p/(1-e^2).
        convertedKeplerianElements( semiMajorAxisIndex ) =
                modifiedEquinoctialElements( semiLatusRectumIndex ) /
                ( getFloatingInteger< ScalarType >( 1 ) - eccentricity * eccentricity );
    }
    // This is (almost) a parabola and the semimajor axis will tend to infinity and be useless.
    else
    {
        // Give semi-latus rectum instead.
        convertedKeplerianElements( semiLatusRectumIndex ) =
                modifiedEquinoctialElements( semiLatusRectumIndex );
    }

    // Compute longitude of ascending node.

    // Compute solution in [-PI,PI] interval with atan2.
    ScalarType longitudeOfAscendingNode = std::atan2( kElement, hElement );

    // Store longitude of ascending node.
    convertedKeplerianElements( longitudeOfAscendingNodeIndex ) = basic_mathematics::
            computeModulo( longitudeOfAscendingNode, getFloatingInteger< ScalarType >( 2 ) * getPi< ScalarType >( ) );

    // Compute inclination.

    // If prograde factor must be positive, otherwise negative.
    ScalarType retrogradeFactor =
            flipSingularityToZeroInclination ? -getFloatingInteger< ScalarType >( 1 ) : getFloatingInteger< ScalarType >( 1 );

    // Calculate magnitude of inclination with retrogradeFactor (derived based on Hintz, 2008).
    ScalarType hSquaredPlusKSquared = std::pow( hElement * hElement + kElement * kElement ,
                                                retrogradeFactor );

    convertedKeplerianElements( inclinationIndex ) =
            getFloatingInteger< ScalarType >( 2 ) * std::atan( std::sqrt( hSquaredPlusKSquared ) );
    // Was: std::atan2(2 * std::sqrt(hSquaredPlusKSquared) , (1.0 - hSquaredPlusKSquared));

    // Compute argument of periapsis.

    // Helper quantity (omega + I * OMEGA ), for argument of periapsis and true anomaly.
    ScalarType argumentOfPeriapsisAndLongitude;

    // If eccentricity is (near) circular
    if ( eccentricity < singularityTolerance )
        // Then the composite argument of periapsis and longitude together cannot be determined
        // because both arguments of the atan2 functions will be zero.
    {
        // Set argument of periapsis to zero.
        convertedKeplerianElements( argumentOfPeriapsisIndex ) = getFloatingInteger< ScalarType >( 0 );

        // Set composite argument to longitude only, since argument of periapsis is now zero.
        argumentOfPeriapsisAndLongitude = retrogradeFactor * longitudeOfAscendingNode;
    }
    else
        // Argument of periapsis can be found from composite argument.
    {
        // Compute composite argument.
        argumentOfPeriapsisAndLongitude = std::atan2( gElement, fElement );

        // Compute argument of periapsis.
        convertedKeplerianElements( argumentOfPeriapsisIndex ) = basic_mathematics::
                computeModulo( argumentOfPeriapsisAndLongitude
                               - retrogradeFactor * longitudeOfAscendingNode,
                               getFloatingInteger< ScalarType >( 2 ) * getPi< ScalarType >( ) );
    }

    // Compute true anomaly.
    convertedKeplerianElements( trueAnomalyIndex ) = basic_mathematics::
            computeModulo( modifiedEquinoctialElements( trueLongitudeIndex )
                           - argumentOfPeriapsisAndLongitude, getFloatingInteger< ScalarType >( 2 ) * getPi< ScalarType >( ) );

    // Return converted elements.
    return convertedKeplerianElements;
}

//! Convert Cartesian to modified equinoctial orbital elements using implicit MEE equation set.
/*!
 * Converts Cartesian to modified equinoctial elements using one of two sets of equations implicitly
 * determined from intermediate inclination. NOTE: This function automatically selects the boolean parameter
 * that defines the location of the singularity, depending on whether the orbit is prograde or retrograde.
 * \param cartesianElements Vector containing Cartesian elements. Order of elements is important!
 *         cartesianElements( 0 ) = x-position coordinate,                                      [m]
 *         cartesianElements( 1 ) = y-position coordinate,                                      [m]
 *         cartesianElements( 2 ) = z-position coordinate,                                      [m]
 *         cartesianElements( 3 ) = x-velocity coordinate,                                    [m/s]
 *         cartesianElements( 4 ) = y-velocity coordinate,                                    [m/s]
 *         cartesianElements( 5 ) = z-velocity coordinate.                                    [m/s]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.
 * \return Converted state in modified equinoctial elements. The order of elements is fixed!
 *         modifiedEquinoctialElements( 0 ) = semi-latus rectum,                                [m]
 *         modifiedEquinoctialElements( 1 ) = f-element,                                        [-]
 *         modifiedEquinoctialElements( 2 ) = g-element,                                        [-]
 *         modifiedEquinoctialElements( 3 ) = h-element,                                        [-]
 *         modifiedEquinoctialElements( 4 ) = k-element,                                        [-]
 *         modifiedEquinoctialElements( 5 ) = true longitude.                                 [rad]
 */
template< typename ScalarType = double >
Eigen::Matrix< ScalarType, 6, 1 > convertCartesianToModifiedEquinoctialElements(
        const Eigen::Matrix< ScalarType, 6, 1 >& cartesianElements,
        const ScalarType centralBodyGravitationalParameter )

{
    // Convert to keplerian elements.
    Eigen::Matrix< ScalarType, 6, 1 > keplerianElements = convertCartesianToKeplerianElements(
                cartesianElements, centralBodyGravitationalParameter );

    // Check whether orbit is retrograde.
    bool flipSingularityToZeroInclination =
            mission_geometry::isOrbitRetrograde( keplerianElements );

    // Convert to modified equinoctial elements.
    return convertKeplerianToModifiedEquinoctialElements(
                keplerianElements, flipSingularityToZeroInclination );
}

template< typename ScalarType = double >
Eigen::Matrix< ScalarType, 6, 1 > convertCartesianToModifiedEquinoctialElementsFromStateFunction(
        const std::function< Eigen::Matrix< ScalarType, 6, 1 >(  ) >& cartesianElementsFunction,
        const std::function< ScalarType( ) > centralBodyGravitationalParameterFunction )
{
    return convertCartesianToModifiedEquinoctialElements(
                cartesianElementsFunction( ), centralBodyGravitationalParameterFunction( ) );
}

//! Convert Cartesian to modified equinoctial orbital elements using explicit MEE equation set.
/*!
 * Converts Cartesian to modified equinoctial elements using one of two sets of equations specified
 * by the user.
 * \param cartesianElements Vector containing Cartesian elements. Order of elements is important!
 *         cartesianElements( 0 ) = x-position coordinate,                                      [m]
 *         cartesianElements( 1 ) = y-position coordinate,                                      [m]
 *         cartesianElements( 2 ) = z-position coordinate,                                      [m]
 *         cartesianElements( 3 ) = x-velocity coordinate,                                    [m/s]
 *         cartesianElements( 4 ) = y-velocity coordinate,                                    [m/s]
 *         cartesianElements( 5 ) = z-velocity coordinate.                                    [m/s]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.
 * \param flipSingularityToZeroInclination Boolean flag to indicate whether the set of equations for
 *          the inclination = 180 degrees (false) or 0 degrees (true) singular case are to be used.
 *          Take note: the same set of equations
 *          is required for conversion back to Cartesian elements to retrieve original state!
 * \return Converted state in modified equinoctial elements. The order of elements is fixed!
 *         modifiedEquinoctialElements( 0 ) = semi-latus rectum,                                [m]
 *         modifiedEquinoctialElements( 1 ) = f-element,                                        [-]
 *         modifiedEquinoctialElements( 2 ) = g-element,                                        [-]
 *         modifiedEquinoctialElements( 3 ) = h-element,                                        [-]
 *         modifiedEquinoctialElements( 4 ) = k-element,                                        [-]
 *         modifiedEquinoctialElements( 5 ) = true longitude.                                 [rad]
 */
template< typename ScalarType = double >
Eigen::Matrix< ScalarType, 6, 1 > convertCartesianToModifiedEquinoctialElements(
        const Eigen::Matrix< ScalarType, 6, 1 >& cartesianElements,
        const ScalarType centralBodyGravitationalParameter,
        const bool flipSingularityToZeroInclination )
{
    // Convert Cartesian to Keplerian to modified equinoctial elements.
    return convertKeplerianToModifiedEquinoctialElements< ScalarType >(
                convertCartesianToKeplerianElements< ScalarType >(
                    cartesianElements, centralBodyGravitationalParameter ),
                flipSingularityToZeroInclination );
}


//! Convert modified equinoctial elements to Cartesian orbital elements.
/*!
 * Converts modified equinoctial elements to Cartesian orbital elements using one of two sets of
 * equations specified by the user. This function first converts to Keplerian elements, and then to Cartesian elements.
 * \param modifiedEquinoctialElements Vector containing modified equinoctial elements. Order of
 * elements is important!
 *         modifiedEquinoctialElements( 0 ) = semi-latus rectum,                                [m]
 *         modifiedEquinoctialElements( 1 ) = f-element,                                        [-]
 *         modifiedEquinoctialElements( 2 ) = g-element,                                        [-]
 *         modifiedEquinoctialElements( 3 ) = h-element,                                        [-]
 *         modifiedEquinoctialElements( 4 ) = k-element,                                        [-]
 *         modifiedEquinoctialElements( 5 ) = true longitude.                                 [rad]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.
 * \param flipSingularityToZeroInclination Boolean flag to indicate whether the set of equations for
 *          the inclination = 180 degrees (false) or 0 degrees (true) singular case are to be used.
 *          Take note: the same set of equations
 *          is required for conversion back to modified equinoctial elements to retrieve original
 *          state!
 * \return Converted state in Cartesian elements. The order of elements is fixed!
 *         cartesianElements( 0 ) = x-position coordinate,                                      [m]
 *         cartesianElements( 1 ) = y-position coordinate,                                      [m]
 *         cartesianElements( 2 ) = z-position coordinate,                                      [m]
 *         cartesianElements( 3 ) = x-velocity coordinate,                                    [m/s]
 *         cartesianElements( 4 ) = y-velocity coordinate,                                    [m/s]
 *         cartesianElements( 5 ) = z-velocity coordinate.                                    [m/s]
 */
template< typename ScalarType = double >
Eigen::Matrix< ScalarType, 6, 1 > convertModifiedEquinoctialToCartesianElementsViaKeplerElements(
        const Eigen::Matrix< ScalarType, 6, 1 >& modifiedEquinoctialElements,
        const ScalarType centralBodyGravitationalParameter,
        const bool flipSingularityToZeroInclination )
// Using unnamed pdf and code archive Bart Rmgens.
{
    using mathematical_constants::getFloatingInteger;

    // Creating output vector.
    Eigen::Matrix< ScalarType, 6, 1 > convertedCartesianElements = Eigen::Matrix< ScalarType, 6, 1 >::Zero( 6 );

    // If the prograde equations are to be used.
    if ( !flipSingularityToZeroInclination )
    {
        // Local storage of elements for ease of access.
        ScalarType semiLatusRectum = modifiedEquinoctialElements( semiLatusRectumIndex );
        ScalarType fElement = modifiedEquinoctialElements( fElementIndex );
        ScalarType gElement = modifiedEquinoctialElements( gElementIndex );
        ScalarType hElement = modifiedEquinoctialElements( hElementIndex );
        ScalarType kElement = modifiedEquinoctialElements( kElementIndex );
        ScalarType trueLongitude = modifiedEquinoctialElements( trueLongitudeIndex );

        // Computing helper parameters.

        // Square-root
        ScalarType squareRootOfGravitationalParameterOverSemiLatusRectum =
                std::sqrt( centralBodyGravitationalParameter / semiLatusRectum );

        // Cosines and sine.
        ScalarType cosineTrueLongitude = std::cos( trueLongitude );
        ScalarType sineTrueLongitude = std::sin( trueLongitude );

        // s, a, w, r from references.
        ScalarType sSquaredParameter = getFloatingInteger< ScalarType >( 1 ) + hElement * hElement + kElement * kElement;
        ScalarType aSquaredParameter = hElement * hElement - kElement*kElement;
        ScalarType wParameter = getFloatingInteger< ScalarType >( 1 ) + fElement * cosineTrueLongitude
                + gElement * sineTrueLongitude;
        ScalarType radius = semiLatusRectum / wParameter;

        // Computing position and storing (using code archive Bart Rmgens and unnamed pdf).
        convertedCartesianElements( xCartesianPositionIndex )
                = radius / sSquaredParameter *
                ( cosineTrueLongitude + aSquaredParameter
                  * cosineTrueLongitude + getFloatingInteger< ScalarType >( 2 ) * hElement * kElement
                  * sineTrueLongitude );
        convertedCartesianElements( yCartesianPositionIndex )
                = radius / sSquaredParameter *
                ( sineTrueLongitude - aSquaredParameter
                  * sineTrueLongitude + getFloatingInteger< ScalarType >( 2 ) * hElement * kElement
                  * cosineTrueLongitude );
        convertedCartesianElements( zCartesianPositionIndex ) = getFloatingInteger< ScalarType >( 2 ) * radius / sSquaredParameter
                * ( hElement * sineTrueLongitude - kElement * cosineTrueLongitude );

        // Computing velocities (these can probably use more optimal computation).
        convertedCartesianElements( xCartesianVelocityIndex ) = -getFloatingInteger< ScalarType >( 1 ) / sSquaredParameter
                * squareRootOfGravitationalParameterOverSemiLatusRectum
                * ( sineTrueLongitude + aSquaredParameter * sineTrueLongitude
                    - getFloatingInteger< ScalarType >( 2 ) * hElement * kElement * cosineTrueLongitude + gElement
                    - getFloatingInteger< ScalarType >( 2 ) * fElement * hElement * kElement + aSquaredParameter * gElement );
        convertedCartesianElements( yCartesianVelocityIndex ) = -getFloatingInteger< ScalarType >( 1 ) / sSquaredParameter
                * squareRootOfGravitationalParameterOverSemiLatusRectum
                * ( -cosineTrueLongitude + aSquaredParameter * cosineTrueLongitude
                    + getFloatingInteger< ScalarType >( 2 ) * hElement * kElement * sineTrueLongitude - fElement
                    + getFloatingInteger< ScalarType >( 2 ) * gElement * hElement * kElement + aSquaredParameter * fElement );
        convertedCartesianElements( zCartesianVelocityIndex ) = getFloatingInteger< ScalarType >( 2 ) / sSquaredParameter
                * squareRootOfGravitationalParameterOverSemiLatusRectum
                * ( hElement * cosineTrueLongitude + kElement * sineTrueLongitude
                    + fElement * hElement + gElement * kElement );
    }
    else // use the indirect transformation via an intermediate Kepler state.
    {
        // Use the Keplerian state as an intermediary step, as there is no direct transformation
        // available in literature (personal derivation forgone for now due to time constraints).
        convertedCartesianElements =
                convertKeplerianToCartesianElements(
                    convertModifiedEquinoctialToKeplerianElements(
                        modifiedEquinoctialElements, flipSingularityToZeroInclination ),
                    centralBodyGravitationalParameter );
    }

    // Return converted set of elements.
    return convertedCartesianElements;
}

//! Convert modified equinoctial elements to Cartesian orbital elements.
/*!
 * Converts modified equinoctial elements to Cartesian orbital elements using one of two sets of
 * equations specified by the user.
 * \param modifiedEquinoctialElements Vector containing modified equinoctial elements. Order of
 * elements is important!
 *         modifiedEquinoctialElements( 0 ) = semi-latus rectum,                                [m]
 *         modifiedEquinoctialElements( 1 ) = f-element,                                        [-]
 *         modifiedEquinoctialElements( 2 ) = g-element,                                        [-]
 *         modifiedEquinoctialElements( 3 ) = h-element,                                        [-]
 *         modifiedEquinoctialElements( 4 ) = k-element,                                        [-]
 *         modifiedEquinoctialElements( 5 ) = true longitude.                                 [rad]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.
 * \param flipSingularityToZeroInclination Boolean flag to indicate whether the set of equations for
 *          the inclination = 180 degrees (false) or 0 degrees (true) singular case are to be used.
 *          Take note: the same set of equations
 *          is required for conversion back to modified equinoctial elements to retrieve original
 *          state!
 * \return Converted state in Cartesian elements. The order of elements is fixed!
 *         cartesianElements( 0 ) = x-position coordinate,                                      [m]
 *         cartesianElements( 1 ) = y-position coordinate,                                      [m]
 *         cartesianElements( 2 ) = z-position coordinate,                                      [m]
 *         cartesianElements( 3 ) = x-velocity coordinate,                                    [m/s]
 *         cartesianElements( 4 ) = y-velocity coordinate,                                    [m/s]
 *         cartesianElements( 5 ) = z-velocity coordinate.                                    [m/s]
 */
template< typename ScalarType = double >
Eigen::Matrix< ScalarType, 6, 1 > convertModifiedEquinoctialToCartesianElements(
        const Eigen::Matrix< ScalarType, 6, 1 >& modifiedEquinoctialElements,
        const ScalarType centralBodyGravitationalParameter,
        const bool flipSingularityToZeroInclination )
{
    if( flipSingularityToZeroInclination )
    {
        throw std::runtime_error(
                    "Error in direct conversion from MEE to Cartesian elements, implementation with singularity at 0 degrees not yet implemented" );
    }
    using mathematical_constants::getFloatingInteger;

    // Retrieve MEE components and compute intermediate quantities
    ScalarType semiLatusRectrum = modifiedEquinoctialElements( semiParameterIndex );
    ScalarType angularMomentumPerUnitGravitationalParameter =
            std::sqrt( semiLatusRectrum / centralBodyGravitationalParameter );

    ScalarType sineTrueLongitude = std::sin( modifiedEquinoctialElements( trueLongitudeIndex ) );
    ScalarType cosineTrueLongitude = std::cos( modifiedEquinoctialElements( trueLongitudeIndex ) );

    ScalarType parameterF = modifiedEquinoctialElements( fElementIndex );
    ScalarType parameterG = modifiedEquinoctialElements( gElementIndex );
    ScalarType parameterH = modifiedEquinoctialElements( hElementIndex );
    ScalarType parameterK = modifiedEquinoctialElements( kElementIndex );

    ScalarType parameterW =
            getFloatingInteger< ScalarType >( 1 ) + parameterF * cosineTrueLongitude + parameterG * sineTrueLongitude;
    ScalarType parameterSSquared = getFloatingInteger< ScalarType >( 1 ) + parameterH * parameterH + parameterK * parameterK;
    ScalarType parameterAlphaSquared = parameterH * parameterH - parameterK * parameterK;

    // Compute Cartesian elements
    Eigen::Matrix< ScalarType, 6, 1 > cartesianElements;
    cartesianElements( 0 ) = cosineTrueLongitude + parameterAlphaSquared * cosineTrueLongitude +
            getFloatingInteger< ScalarType >( 2 ) * parameterH * parameterK * sineTrueLongitude;
    cartesianElements( 1 ) = ( sineTrueLongitude - parameterAlphaSquared * sineTrueLongitude +
              getFloatingInteger< ScalarType >( 2 ) * parameterH * parameterK * cosineTrueLongitude );
    cartesianElements( 2 ) =
            getFloatingInteger< ScalarType >( 2 ) * ( parameterH * sineTrueLongitude - parameterK * cosineTrueLongitude );
    cartesianElements.segment( 0, 3 ) *= semiLatusRectrum / ( parameterW * parameterSSquared );

    cartesianElements( 3 ) =
            -( sineTrueLongitude + parameterAlphaSquared * sineTrueLongitude -
               getFloatingInteger< ScalarType >( 2 ) * parameterH * parameterK * cosineTrueLongitude + parameterG -
               getFloatingInteger< ScalarType >( 2 ) * parameterF * parameterH * parameterK +
               parameterAlphaSquared * parameterG );
    cartesianElements( 4 ) = -( -cosineTrueLongitude + parameterAlphaSquared * cosineTrueLongitude +
               getFloatingInteger< ScalarType >( 2 ) * parameterH * parameterK * sineTrueLongitude -
               parameterF + getFloatingInteger< ScalarType >( 2 ) * parameterG * parameterH * parameterK +
               parameterAlphaSquared * parameterF );
    cartesianElements( 5 ) =
            getFloatingInteger< ScalarType >( 2 ) * ( parameterH * cosineTrueLongitude + parameterK * sineTrueLongitude +
                                                      parameterF * parameterH +  parameterG * parameterK );
    cartesianElements.segment( 3, 3 ) *=
            getFloatingInteger< ScalarType >( 1 ) / ( parameterSSquared * angularMomentumPerUnitGravitationalParameter );

    return cartesianElements;


}

} // namespace orbital_element_conversions

} // namespace tudat

#endif // TUDAT_MODIFIED_EQUINOCTIAL_ELEMENT_CONVERSIONS_H

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
 *      Verified Interval Propagation, Bart Rmgens; Delft (2011). Code archive.
 *      Modified Equinoctial Orbital Elements, author unknown;
 *          http://www.cdeagle.com/pdf/mee.pdf (2010?).
 *      Survey of Orbital Element Sets, Gerald R. Hintz; Journal of Guidance, Control and
 *          Dynamics (2008, Vol. 31 - Nr. 3).
 *      Code archive, E. Heeren (fellow Tudat developer).
 *
 */

#include <cmath>


#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/missionGeometry.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/modifiedEquinoctialElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"

namespace tudat
{

namespace orbital_element_conversions
{


//! Convert Keplerian to modified equinoctial orbital elements using implicit MEE equation set.
Eigen::Vector6d convertKeplerianToModifiedEquinoctialElements(
        const Eigen::Vector6d& keplerianElements )
// Based on Hintz, 2008.
{
    // Check if orbit is retrograde
    bool avoidSingularityAtPiInclination =
            mission_geometry::isOrbitRetrograde( keplerianElements );

    // Convert to modified equinoctial elements.
    return convertKeplerianToModifiedEquinoctialElements( keplerianElements,
                                                          avoidSingularityAtPiInclination );
}

//! Convert Keplerian to modified equinoctial orbital elements using MEE explicit equation set.
Eigen::Vector6d convertKeplerianToModifiedEquinoctialElements(
        const Eigen::Vector6d& keplerianElements,
        const bool avoidSingularityAtPiInclination )
// Based on Hintz, 2008.
{
    using mathematical_constants::PI;

    // Declaring eventual output vector.
    Eigen::Vector6d modifiedEquinoctialState( 6 );

    // Compute semi-latus rectum.
    double singularityTolerance = 1.0e-15; // Based on tolerance chosen in
                                           // orbitalElementConversions.cpp in Tudat Core.

    // Extracting eccentricity for ease of referencing.
    double eccentricity = keplerianElements( eccentricityIndex );

    // If e is (very near) one, then semi-major axis is undefined and thus the first kepler
    // element is the semi-latus rectum.
    if ( std::fabs( eccentricity - 1.0 ) < singularityTolerance )
    {
        modifiedEquinoctialState( semiLatusRectumIndex )
                = keplerianElements( semiLatusRectumIndex );
    }
    else // eccentricity is significantly away from singularity and can be computed with p=a(1-e).
    {
        modifiedEquinoctialState( semiLatusRectumIndex ) = keplerianElements( semiMajorAxisIndex )
                                                           * ( 1.0 - eccentricity * eccentricity );
    }

    // Determine prograde-ness, as other five parameters depend on that fact.
    double inclination = keplerianElements( inclinationIndex );

    // If inclination is outside range [0,PI].
    if ( ( inclination < 0.0 ) || ( inclination > PI ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "Inclination is expected in range [0," << PI << "]\n"
                     << "Specified inclination: " << inclination << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( errorMessage.str( ) );
    }
    //Else, nothing wrong and continue.

    // Compute set dependant helper parameters.
    double argumentOfPeriapsisAndAscendingNode = 0.0;
    double tangentOfHalfInclination = 0.0;

    // If normal set of equations is to be used (i.e. not inverse set for singularity).
    if ( !avoidSingularityAtPiInclination )
    {
        // Add (+) argument of periapsis and longitude of ascending node.
        argumentOfPeriapsisAndAscendingNode = keplerianElements( argumentOfPeriapsisIndex )
                + keplerianElements( longitudeOfAscendingNodeIndex );

        // Take the tangent of half the inclination.
        tangentOfHalfInclination = std::tan( inclination / 2.0 );
    }
    else // the orbit is retrograde.
    {
        // Subtract (-) longitude of ascending node from argument of periapsis.
        argumentOfPeriapsisAndAscendingNode = keplerianElements( argumentOfPeriapsisIndex )
                - keplerianElements( longitudeOfAscendingNodeIndex );

        // Take the inverse of the tangent of half the inclination to avoid singularity at tan
        // PI/2.
        tangentOfHalfInclination = 1.0 / std::tan( inclination / 2.0 );
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
                + keplerianElements( trueAnomalyIndex ), 2.0 * PI );

    // Give back result.
    return modifiedEquinoctialState;
}

//! Convert modified equinoctial to Keplerian orbital elements.
Eigen::Vector6d convertModifiedEquinoctialToKeplerianElements(
        const Eigen::Vector6d& modifiedEquinoctialElements,
        const bool avoidSingularityAtPiInclination )
// Using unknown source pdf, code archive E. Heeren and personal derivation based on Hintz 2008.
{
    using mathematical_constants::PI;

    // Declaration of output vector.
    Eigen::Vector6d convertedKeplerianElements = Eigen::Vector6d::Zero( 6 );

    // for ease of referencing, almost all modified equinoctial elements.
    double fElement = modifiedEquinoctialElements( fElementIndex );
    double gElement = modifiedEquinoctialElements( gElementIndex );
    double hElement = modifiedEquinoctialElements( hElementIndex );
    double kElement = modifiedEquinoctialElements( kElementIndex );

    // Tolerance for singularities.
    double singularityTolerance = 1.0e-15;

    // Compute eccentricity.
    double eccentricity = std::sqrt( fElement * fElement + gElement * gElement );
    convertedKeplerianElements( eccentricityIndex ) = eccentricity;

    // Compute semi-major axis.
    // If eccentricity is not near-parabolic.
    if ( std::fabs( eccentricity - 1.0 ) > singularityTolerance )
    {
        // Use semi-latus rectum and eccentricity to calculate semi-major axis with a=p/(1-e^2).
        convertedKeplerianElements( semiMajorAxisIndex ) =
                modifiedEquinoctialElements( semiLatusRectumIndex ) /
                ( 1.0 - eccentricity * eccentricity );
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
    double longitudeOfAscendingNode = std::atan2( kElement, hElement );

    // Store longitude of ascending node.
    convertedKeplerianElements( longitudeOfAscendingNodeIndex ) = basic_mathematics::
            computeModulo( longitudeOfAscendingNode, 2.0 * PI );

    // Compute inclination.

    // If prograde factor must be positive, otherwise negative.
    double retrogradeFactor = avoidSingularityAtPiInclination ? -1.0 : 1.0;

    // Calculate magnitude of inclination with retrogradeFactor (derived based on Hintz, 2008).
    double hSquaredPlusKSquared = std::pow( hElement * hElement + kElement * kElement ,
                                            retrogradeFactor );

    convertedKeplerianElements( inclinationIndex ) =
            2.0 * std::atan( std::sqrt( hSquaredPlusKSquared ) );
    // Was: std::atan2(2.0 * std::sqrt(hSquaredPlusKSquared) , (1.0 - hSquaredPlusKSquared));

    // Compute argument of periapsis.

    // Helper quantity (omega + I * OMEGA ), for argument of periapsis and true anomaly.
    double argumentOfPeriapsisAndLongitude;

    // If eccentricity is (near) circular
    if ( eccentricity < singularityTolerance )
        // Then the composite argument of periapsis and longitude together cannot be determined
        // because both arguments of the atan2 functions will be zero.
    {
        // Set argument of periapsis to zero.
        convertedKeplerianElements( argumentOfPeriapsisIndex ) = 0.0;

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
                               - retrogradeFactor * longitudeOfAscendingNode, 2.0 * PI );
    }

    // Compute true anomaly.
    convertedKeplerianElements( trueAnomalyIndex ) = basic_mathematics::
            computeModulo( modifiedEquinoctialElements( trueLongitudeIndex )
                           - argumentOfPeriapsisAndLongitude, 2.0 * PI );

    // Return converted elements.
    return convertedKeplerianElements;
}

//! Convert Cartesian to modified equinoctial orbital elements using implicit MEE equation set.
Eigen::Vector6d convertCartesianToModifiedEquinoctialElements(
        const Eigen::Vector6d& cartesianElements,
        const double centralBodyGravitationalParameter )
{
    using mathematical_constants::PI;

    // Convert to keplerian elements.
    Eigen::Vector6d keplerianElements = convertCartesianToKeplerianElements(
                cartesianElements, centralBodyGravitationalParameter );

    // Check whether orbit is retrograde.
    bool avoidSingularityAtPiInclination =
            mission_geometry::isOrbitRetrograde( keplerianElements );

    // Convert to modified equinoctial elements.
    return convertKeplerianToModifiedEquinoctialElements(
                keplerianElements, avoidSingularityAtPiInclination );
}

//! Convert Cartesian to modified equinoctial orbital elements using explicit MEE equation set.
Eigen::Vector6d convertCartesianToModifiedEquinoctialElements(
        const Eigen::Vector6d& cartesianElements,
        const double centralBodyGravitationalParameter,
        const bool avoidSingularityAtPiInclination )
{
    // Convert Cartesian to Keplerian to modified equinoctial elements.
    return convertKeplerianToModifiedEquinoctialElements(
                convertCartesianToKeplerianElements(
                    cartesianElements, centralBodyGravitationalParameter ),
                avoidSingularityAtPiInclination );
}

//! Convert Modified Equinoctial Elements to Cartesian Elements.
Eigen::Vector6d convertModifiedEquinoctialToCartesianElements(
        const Eigen::Vector6d& modifiedEquinoctialElements,
        const double centralBodyGravitationalParameter,
        const bool avoidSingularityAtPiInclination )
// Using unnamed pdf and code archive Bart Rmgens.
{
    // Creating output vector.
    Eigen::Vector6d convertedCartesianElements = Eigen::Vector6d::Zero( 6 );

    // If the prograde equations are to be used.
    if ( !avoidSingularityAtPiInclination )
    {
        // Local storage of elements for ease of access.
        double semiLatusRectum = modifiedEquinoctialElements( semiLatusRectumIndex );
        double fElement = modifiedEquinoctialElements( fElementIndex );
        double gElement = modifiedEquinoctialElements( gElementIndex );
        double hElement = modifiedEquinoctialElements( hElementIndex );
        double kElement = modifiedEquinoctialElements( kElementIndex );
        double trueLongitude = modifiedEquinoctialElements( trueLongitudeIndex );

        // Computing helper parameters.

        // Square-root
        double squareRootOfGravitationalParameterOverSemiLatusRectum =
                std::sqrt( centralBodyGravitationalParameter / semiLatusRectum );

        // Cosines and sine.
        double cosineTrueLongitude = std::cos( trueLongitude );
        double sineTrueLongitude = std::sin( trueLongitude );

        // s, a, w, r from references.
        double sSquaredParameter = 1.0 + hElement * hElement + kElement * kElement;
        double aSquaredParameter = hElement * hElement - kElement*kElement;
        double wParameter = 1.0 + fElement * cosineTrueLongitude
                + gElement * sineTrueLongitude;
        double radius = semiLatusRectum / wParameter;

        // Computing position and storing (using code archive Bart Rmgens and unnamed pdf).
        convertedCartesianElements( xCartesianPositionIndex )
                = radius / sSquaredParameter * ( cosineTrueLongitude + aSquaredParameter
                                                 * cosineTrueLongitude + 2.0 * hElement * kElement
                                                 * sineTrueLongitude );
        convertedCartesianElements( yCartesianPositionIndex )
                = radius / sSquaredParameter * ( sineTrueLongitude - aSquaredParameter
                                                 * sineTrueLongitude + 2.0 * hElement * kElement
                                                 * cosineTrueLongitude );
        convertedCartesianElements( zCartesianPositionIndex ) = 2.0 * radius / sSquaredParameter
                * ( hElement * sineTrueLongitude - kElement * cosineTrueLongitude );

        // Computing velocities (these can probably use more optimal computation).
        convertedCartesianElements( xCartesianVelocityIndex ) = -1.0 / sSquaredParameter
                * squareRootOfGravitationalParameterOverSemiLatusRectum
                * ( sineTrueLongitude + aSquaredParameter * sineTrueLongitude
                    - 2.0 * hElement * kElement * cosineTrueLongitude + gElement
                    - 2.0 * fElement * hElement * kElement + aSquaredParameter * gElement );
        convertedCartesianElements( yCartesianVelocityIndex ) = -1.0 / sSquaredParameter
                * squareRootOfGravitationalParameterOverSemiLatusRectum
                * ( -cosineTrueLongitude + aSquaredParameter * cosineTrueLongitude
                    + 2.0 * hElement * kElement * sineTrueLongitude - fElement
                    + 2.0 * gElement * hElement * kElement + aSquaredParameter * fElement );
        convertedCartesianElements( zCartesianVelocityIndex ) = 2.0 / sSquaredParameter
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
                        modifiedEquinoctialElements, avoidSingularityAtPiInclination ),
                    centralBodyGravitationalParameter );
    }

    // Return converted set of elements.
    return convertedCartesianElements;
}

} // namespace orbital_element_conversions

} // namespace tudat

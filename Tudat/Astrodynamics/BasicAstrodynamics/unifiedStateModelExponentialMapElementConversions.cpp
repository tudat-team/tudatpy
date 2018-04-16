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
 *      Vittaldev, V. (2010). The unified state model: Derivation and application in astrodynamics
 *          and navigation. Master thesis, Delft University of Technology.
 *      Facchinelli, M. (2018). Aerobraking Guidance, Navigation and Control.
 *          Master thesis, Delft University of Technology.
 */

#include <cmath>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/missionGeometry.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unifiedStateModelExponentialMapElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"

namespace tudat
{

namespace orbital_element_conversions
{

//! Convert Keplerian elements to unified state model elements with exponential map.
Eigen::Vector6d convertKeplerianToUnifiedStateModelExponentialMapElements(
        const Eigen::Vector6d& keplerianElements,
        const double centralBodyGravitationalParameter )
{
    using mathematical_constants::PI;

    // Declaring eventual output vector.
    Eigen::Vector6d convertedUnifiedStateModelElements = Eigen::Vector6d::Zero( );

    // Define the tolerance of a singularity
    double singularityTolerance = 20.0 * std::numeric_limits< double >::epsilon( );

    // If eccentricity is outside range [0,inf)
    if ( keplerianElements( eccentricityIndex ) < 0.0 )
    {
        //Define the error message
        std::stringstream errorMessage;
        errorMessage << "Eccentricity is expected in range [0,inf)\n"
                     << "Specified eccentricity: " << keplerianElements( eccentricityIndex ) << std::endl;

        // Throw exception
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }

    // If inclination is outside range [0,PI]
    if ( ( keplerianElements( inclinationIndex ) < 0.0 ) || ( keplerianElements( inclinationIndex ) > PI ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "Inclination is expected in range [0," << PI << "]\n"
                     << "Specified inclination: " << keplerianElements( inclinationIndex ) << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }

    // If argument of pericenter is outside range [0,2.0 * PI]
    if ( ( keplerianElements( argumentOfPeriapsisIndex ) < 0.0 ) || ( keplerianElements( argumentOfPeriapsisIndex ) >
                                                                      2.0 * PI ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "Argument of periapsis is expected in range [0," << 2.0 * PI << "]\n"
                     << "Specified inclination: " << keplerianElements( argumentOfPeriapsisIndex ) << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }

    // If right ascension of ascending node is outside range [0,2.0 * PI]
    if ( ( keplerianElements( longitudeOfAscendingNodeIndex ) < 0.0 ) ||
         ( keplerianElements( longitudeOfAscendingNodeIndex ) > 2.0 * PI ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "RAAN is expected in range [0," << 2.0 * PI << "]\n"
                     << "Specified inclination: " << keplerianElements( longitudeOfAscendingNodeIndex ) << " rad."
                     << std::endl;

        // Throw exception.
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }

    // If true anomaly is outside range [0,2.0 * PI]
    if ( ( keplerianElements( trueAnomalyIndex ) < 0.0 ) || ( keplerianElements( trueAnomalyIndex ) > 2.0 * PI ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "True anomaly is expected in range [0," << 2.0 * PI << "]\n"
                     << "Specified inclination: " << keplerianElements( trueAnomalyIndex ) << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }

    // If inclination is zero and the right ascension of ascending node is non-zero
    if ( ( std::fabs( keplerianElements( inclinationIndex ) ) < singularityTolerance ) &&
         ( std::fabs( keplerianElements( longitudeOfAscendingNodeIndex ) ) > singularityTolerance ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "When the inclination is zero, the right ascending node should be zero by definition\n"
                     << "Specified right ascension of ascending node: " <<
                        keplerianElements( longitudeOfAscendingNodeIndex ) << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }

    // If eccentricity is zero and the argument of pericenter is non-zero
    if ( ( std::fabs( keplerianElements( eccentricityIndex ) ) < singularityTolerance ) &&
         ( std::fabs( keplerianElements( argumentOfPeriapsisIndex ) ) > singularityTolerance ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "When the eccentricity is zero, the argument of pericenter should be zero by definition\n"
                     << "Specified argument of pericenter: " <<
                        keplerianElements( argumentOfPeriapsisIndex ) << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }

    // If semi-major axis is negative and the eccentricity is smaller or equal to one
    if ( ( keplerianElements( semiMajorAxisIndex ) < 0.0 ) && ( keplerianElements( eccentricityIndex ) <= 1.0 ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "When the semi-major axis is negative, the eccentricity should be larger than one\n"
                     << "Specified semi-major axis: " << keplerianElements( semiMajorAxisIndex ) << " m.\n"
                     << "Specified eccentricity: " << keplerianElements( eccentricityIndex ) << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }

    // If semi-major axis is positive and the eccentricity is larger than one
    if ( ( keplerianElements( semiMajorAxisIndex ) > 0.0 ) && ( keplerianElements( eccentricityIndex ) > 1.0 ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "When the semi-major axis is positive, the eccentricity should be smaller than or equal to one\n"
                     << "Specified semi-major axis: " << keplerianElements( semiMajorAxisIndex ) << " m.\n"
                     << "Specified eccentricity: " << keplerianElements( eccentricityIndex ) << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }
    //Else, nothing wrong and continue

    // Compute the C hodograph element of the unified state model
    if ( std::fabs( keplerianElements( eccentricityIndex ) - 1.0) < singularityTolerance )
        // parabolic orbit -> semi-major axis is not defined
    {
        convertedUnifiedStateModelElements( CHodographExponentialMapIndex ) =
                std::sqrt( centralBodyGravitationalParameter / keplerianElements( semiLatusRectumIndex ) );
    }
    else
    {
        convertedUnifiedStateModelElements( CHodographExponentialMapIndex ) =
                std::sqrt( centralBodyGravitationalParameter / ( keplerianElements( semiMajorAxisIndex )
                                                                 * ( 1 - keplerianElements( eccentricityIndex ) *
                                                                     keplerianElements( eccentricityIndex ) ) ) );
    }

    // Calculate the additional R hodograph parameter
    double RHodographElement = keplerianElements( eccentricityIndex ) *
            convertedUnifiedStateModelElements( CHodographExponentialMapIndex );

    // Compute the Rf1 hodograph element of the unified state model
    convertedUnifiedStateModelElements( Rf1HodographExponentialMapIndex ) =
            - RHodographElement * std::sin( keplerianElements( longitudeOfAscendingNodeIndex )
                                            + keplerianElements( argumentOfPeriapsisIndex ) );

    // Compute the Rf2 hodograph element of the unified state model
    convertedUnifiedStateModelElements( Rf2HodographExponentialMapIndex ) =
            RHodographElement * std::cos( keplerianElements( longitudeOfAscendingNodeIndex )
                                          + keplerianElements( argumentOfPeriapsisIndex ) );

    // Calculate the additional elements
    double argumentOfLatitude = keplerianElements( argumentOfPeriapsisIndex ) +
            keplerianElements( trueAnomalyIndex ); // also called u
    double rightAscensionOfLatitude = argumentOfLatitude +
            keplerianElements( longitudeOfAscendingNodeIndex ); // also called lambda

    // Compute magnitude of exponential map
    double arccosineArgument = std::cos( 0.5 * keplerianElements( inclinationIndex ) ) *
            std::cos( 0.5 * rightAscensionOfLatitude );
    if ( ( std::fabs( arccosineArgument ) - 1.0 ) > 0.0 )
    {
        // Make sure that the cosine does not exceed 1.0 in magnitude
        arccosineArgument = ( arccosineArgument > 0.0 ) ? 1.0 : - 1.0;
    }
    double exponentialMapMagnitude = 2.0 * std::acos( arccosineArgument );

    // Check for singularity
    if ( std::fabs( exponentialMapMagnitude ) < singularityTolerance ||
         std::fabs( exponentialMapMagnitude - 2.0 * PI ) < singularityTolerance )
    {
        // If rotation angle is zero (or 2 PI), the exponential map vector is the zero vector
        convertedUnifiedStateModelElements.segment( e1ExponentialMapIndex, 3 ) = Eigen::Vector3d::Zero( );
    }
    else
    {
        // Find the common multiplication factor to the vector elements
        double multiplicationFactor = exponentialMapMagnitude / std::sin( 0.5 * exponentialMapMagnitude );

        // Compute the e1 exponential map of the unified state model
        convertedUnifiedStateModelElements( e1ExponentialMapIndex ) =
                multiplicationFactor * std::sin( 0.5 * keplerianElements( inclinationIndex ) ) *
                std::cos( 0.5 * ( keplerianElements( longitudeOfAscendingNodeIndex ) - argumentOfLatitude ) );

        // Compute the e2 exponential map of the unified state model
        convertedUnifiedStateModelElements( e2ExponentialMapIndex ) =
                multiplicationFactor * std::sin( 0.5 * keplerianElements( inclinationIndex ) ) *
                std::sin( 0.5 * ( keplerianElements( longitudeOfAscendingNodeIndex ) - argumentOfLatitude ) );

        // Compute the e3 exponential map of the unified state model
        convertedUnifiedStateModelElements( e3ExponentialMapIndex ) =
                multiplicationFactor * std::cos( 0.5 * keplerianElements( inclinationIndex ) ) *
                std::sin( 0.5 * rightAscensionOfLatitude );
    }

    // Give back result
    return convertedUnifiedStateModelElements;
}

//! Convert unified state model elements with exponential map to Keplerian elements.
Eigen::Vector6d convertUnifiedStateModelExponentialMapToKeplerianElements(
        const Eigen::Vector6d& unifiedStateModelElements,
        const double centralBodyGravitationalParameter )
{
    using mathematical_constants::PI;

    // Declaring eventual output vector.
    Eigen::Vector6d convertedKeplerianElements = Eigen::Vector6d::Zero( );

    // Define the tolerance of a singularity
    double singularityTolerance = 20.0 * std::numeric_limits< double >::epsilon( );

    // Compute auxiliary parameters
    Eigen::Vector3d exponentialMapVector = unifiedStateModelElements.segment( e1ExponentialMapIndex, 3 );
    double exponentialMapMagnitude = exponentialMapVector.norm( ); // magnitude of exponential map, also called xi
    Eigen::Vector3d eulerEigenaxisVector = Eigen::Vector3d::Zero( ); // Euler eigenaxis vector (unit vector)

    // Compute right ascension of latitude
    double rightAscensionOfLatitude = 0.0;
    if ( exponentialMapMagnitude < singularityTolerance )
    {
        // When exponential map is zero, all Keplerian angles are also zero
        convertedKeplerianElements( inclinationIndex ) = 0.0;
        convertedKeplerianElements( argumentOfPeriapsisIndex ) = 0.0;
        convertedKeplerianElements( longitudeOfAscendingNodeIndex ) = 0.0;
        convertedKeplerianElements( trueAnomalyIndex ) = 0.0;

        // Naturally, lambda is also zero
        rightAscensionOfLatitude = 0.0;
    }
    else
    {
        // Find Euler eigenaxis vector
        eulerEigenaxisVector = exponentialMapVector / exponentialMapMagnitude;

        // Find lambda
        rightAscensionOfLatitude = 2.0 * std::atan2( eulerEigenaxisVector( e3ExponentialMapIndex - 3 ) *
                                                     std::sin( 0.5 * exponentialMapMagnitude ),
                                                     std::cos( 0.5 * exponentialMapMagnitude ) );
    }

    // Trigonometric values of right ascension of latitude
    double cosineLambda = std::cos( rightAscensionOfLatitude );
    double sineLambda = std::sin( rightAscensionOfLatitude );

    // Compute auxiliary parameters auxiliaryParameter1 and auxiliaryParameter2
    double auxiliaryParameter1 = unifiedStateModelElements( Rf1HodographExponentialMapIndex ) * cosineLambda +
            unifiedStateModelElements( Rf2HodographExponentialMapIndex ) * sineLambda;
    double auxiliaryParameter2 = unifiedStateModelElements( CHodographExponentialMapIndex ) -
            unifiedStateModelElements( Rf1HodographExponentialMapIndex ) * sineLambda +
            unifiedStateModelElements( Rf2HodographExponentialMapIndex ) * cosineLambda;

    // Compute auxiliary R hodograph parameter
    double RHodographElement = std::sqrt( unifiedStateModelElements( Rf1HodographExponentialMapIndex ) *
                                          unifiedStateModelElements( Rf1HodographExponentialMapIndex ) +
                                          unifiedStateModelElements( Rf2HodographExponentialMapIndex ) *
                                          unifiedStateModelElements( Rf2HodographExponentialMapIndex ) );

    // Compute eccentricity
    convertedKeplerianElements( eccentricityIndex ) =
            RHodographElement / unifiedStateModelElements( CHodographExponentialMapIndex );

    // Compute semi-major axis or, in case of a parabolic orbit, the semi-latus rectum.
    if ( std::fabs( convertedKeplerianElements( eccentricityIndex ) - 1.0 ) < singularityTolerance )
        // parabolic orbit -> semi-major axis is not defined. Use semi-latus rectum instead.
    {
        convertedKeplerianElements( semiLatusRectumIndex ) = centralBodyGravitationalParameter /
                ( unifiedStateModelElements( CHodographExponentialMapIndex ) *
                  unifiedStateModelElements( CHodographExponentialMapIndex ) );
    }
    else
    {
        convertedKeplerianElements( semiMajorAxisIndex ) =
                centralBodyGravitationalParameter /
                ( std::pow( unifiedStateModelElements( CHodographExponentialMapIndex ), 2 ) *
                  ( 1 - std::pow( convertedKeplerianElements( eccentricityIndex ), 2 ) ) );
    }

    // Continue only if angles have not been computed yet
    if ( exponentialMapMagnitude < singularityTolerance )
    {
        // Give back result
        return convertedKeplerianElements;
    }
    else
    {
        // Compute inclination
        convertedKeplerianElements( inclinationIndex ) = std::acos(
                    ( std::pow( eulerEigenaxisVector( e3ExponentialMapIndex - 3 ), 2 ) - 1.0 ) *
                    ( 1.0 - std::cos( exponentialMapMagnitude ) ) + 1.0 );
        // this acos is always defined correctly because the inclination is always below pi rad

        // Find sine and cosine of longitude of ascending node separately
        double sineOmega = eulerEigenaxisVector( e1ExponentialMapIndex - 3 ) * eulerEigenaxisVector( e3ExponentialMapIndex - 3 ) *
                ( 1 - std::cos( exponentialMapMagnitude ) ) + eulerEigenaxisVector( e2ExponentialMapIndex - 3 ) *
                std::sin( exponentialMapMagnitude );
        double cosineOmega = - eulerEigenaxisVector( e2ExponentialMapIndex - 3 ) * eulerEigenaxisVector( e3ExponentialMapIndex - 3 ) *
                ( 1 - std::cos( exponentialMapMagnitude ) ) + eulerEigenaxisVector( e1ExponentialMapIndex - 3 ) *
                std::sin( exponentialMapMagnitude );
        double denominator = std::sqrt( cosineOmega * cosineOmega +
                                        sineOmega * sineOmega );

        // Compute longitude of ascending node
        if ( std::fabs( std::fabs( convertedKeplerianElements( inclinationIndex ) ) - PI ) < singularityTolerance )
            // pure-prograde or pure-retrograde orbit
        {
            // Define the error message
            std::stringstream errorMessage;
            errorMessage << "Pure-retrograde orbit (inclination = pi).\n"
                         << "Unified state model elements cannot be transformed to Kepler elements." << std::endl;

            // Throw exception
            throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
        }
        else if ( std::fabs( denominator ) < singularityTolerance )
            // null denominator, find work-around
        {
            convertedKeplerianElements( longitudeOfAscendingNodeIndex ) = 0.0; // by definition
        }
        else
        {
            // Compute longitude of ascending node
            convertedKeplerianElements( longitudeOfAscendingNodeIndex ) = std::atan2(
                        sineOmega / denominator, cosineOmega / denominator );

            // Round off small values of the right ascension of ascending node to zero
            if ( std::fabs( convertedKeplerianElements( longitudeOfAscendingNodeIndex ) ) < singularityTolerance )
            {
                convertedKeplerianElements( longitudeOfAscendingNodeIndex ) = 0.0;
            }

            // Ensure the longitude of ascending node is positive
            while ( convertedKeplerianElements( longitudeOfAscendingNodeIndex ) < 0.0 )
                // Because of the previous if statement, if the longitude of ascending node is smaller than 0, it will
                // always be smaller than -singularityTolerance
            {
                convertedKeplerianElements( longitudeOfAscendingNodeIndex ) += 2.0 * PI;
            }
        }

        // Compute true anomaly and argument of periapsis
        if ( std::fabs( RHodographElement ) < singularityTolerance ) // circular orbit
        {
            convertedKeplerianElements( argumentOfPeriapsisIndex ) = 0.0; // by definition
            convertedKeplerianElements( trueAnomalyIndex ) =
                    rightAscensionOfLatitude - convertedKeplerianElements( longitudeOfAscendingNodeIndex );

            // Round off small theta to zero
            if ( std::fabs( convertedKeplerianElements( trueAnomalyIndex ) ) < singularityTolerance )
            {
                convertedKeplerianElements( trueAnomalyIndex ) = 0.0;
            }

            // Ensure the true anomaly is positive
            while ( convertedKeplerianElements( trueAnomalyIndex ) < 0.0 )
                // Because of the previous if statement, if the true anomaly is smaller than zero, it will always be smaller than
                // -singularityTolerance
            {
                convertedKeplerianElements( trueAnomalyIndex ) += 2.0 * PI;
            }
        }
        else
        {
            convertedKeplerianElements( trueAnomalyIndex ) =
                    std::atan2( ( auxiliaryParameter1 / RHodographElement ),
                                ( ( auxiliaryParameter2 - unifiedStateModelElements( CHodographExponentialMapIndex ) )
                                  / RHodographElement ) );

            // Round off small theta to zero
            if ( std::fabs( convertedKeplerianElements( trueAnomalyIndex ) ) < singularityTolerance )
            {
                convertedKeplerianElements( trueAnomalyIndex ) = 0.0;
            }

            // Ensure the true anomaly is positive
            while ( convertedKeplerianElements( trueAnomalyIndex ) < 0.0 )
                // Because of the previous if statement, if the true anomaly is smaller than zero, it will always
                // be smaller than -singularityTolerance
            {
                convertedKeplerianElements( trueAnomalyIndex ) += 2.0 * PI;
            }

            convertedKeplerianElements( argumentOfPeriapsisIndex ) =
                    rightAscensionOfLatitude -
                    convertedKeplerianElements( longitudeOfAscendingNodeIndex ) -
                    convertedKeplerianElements( trueAnomalyIndex );

            // Round off small omega to zero
            if ( std::fabs( convertedKeplerianElements( argumentOfPeriapsisIndex ) ) < singularityTolerance )
            {
                convertedKeplerianElements( argumentOfPeriapsisIndex ) = 0.0;
            }

            // Ensure the argument of periapsis is positive
            while ( convertedKeplerianElements( argumentOfPeriapsisIndex ) < 0.0 )
                // Because of the previous if statement, if the argument of pericenter is smaller than zero,
                // it will be smaller than -singularityTolerance
            {
                convertedKeplerianElements( argumentOfPeriapsisIndex ) += 2.0 * PI;
            }
        }

        // Give back result
        return convertedKeplerianElements;
    }
}

//! Convert Cartesian elements to unified state model elements with exponential map.
Eigen::Vector6d convertCartesianToUnifiedStateModelExponentialMapElements(
        const Eigen::Vector6d& cartesianElements,
        const double centralBodyGravitationalParameter )
{
    using mathematical_constants::PI;

    // Declaring eventual output vector.
    Eigen::Vector6d convertedUnifiedStateModelExponentialMapElements = Eigen::Vector6d::Zero( );

    // Define the tolerance of a singularity
    double singularityTolerance = 20.0 * std::numeric_limits< double >::epsilon( );

    // Convert Cartesian to USM7
    Eigen::Vector7d unifiedStateModelQuaternionsElements =
            convertCartesianToUnifiedStateModelQuaternionsElements( cartesianElements,
                                                                    centralBodyGravitationalParameter );

    // Extract quaternions
    Eigen::Vector3d epsilonQuaternionVector =
            unifiedStateModelQuaternionsElements.segment( epsilon1QuaternionIndex, 3 );
    double etaQuaternionParameter = unifiedStateModelQuaternionsElements( etaQuaternionIndex );

    // Convert quaternions to exponential map (or SEM)
    double exponentialMapMagnitude = 2.0 * std::acos( etaQuaternionParameter );
    if ( std::fabs( exponentialMapMagnitude ) < singularityTolerance ||
         std::fabs( exponentialMapMagnitude - 2.0 * PI ) < singularityTolerance )
    {
        // If rotation angle is zero (or 2 PI), the exponential map vector is the zero vector
        convertedUnifiedStateModelExponentialMapElements.segment( e1ExponentialMapIndex, 3 ) =
                Eigen::Vector3d::Zero( );
    }
    else
    {
        // Find Euler eigenaxis vector
        Eigen::Vector3d eulerEigenaxisVector = epsilonQuaternionVector /
                std::sin( 0.5 * exponentialMapMagnitude );
        convertedUnifiedStateModelExponentialMapElements.segment( e1ExponentialMapIndex, 3 ) =
                ( exponentialMapMagnitude > PI ) ?
                    - ( 2.0 * PI - exponentialMapMagnitude ) * eulerEigenaxisVector : // shadow exponential map
                    exponentialMapMagnitude * eulerEigenaxisVector; // exponential map
    }

    // Add other elements to USMEM vector
    convertedUnifiedStateModelExponentialMapElements.segment( CHodographExponentialMapIndex, 3 ) =
            unifiedStateModelQuaternionsElements.segment( CHodographQuaternionIndex, 3 );

    // Give back result
    return convertedUnifiedStateModelExponentialMapElements;
}

//! Convert unified state model elements with exponential map to Cartesian elements.
Eigen::Vector6d convertUnifiedStateModelExponentialMapToCartesianElements(
        const Eigen::Vector6d& unifiedStateModelExponentialMapElements,
        const double centralBodyGravitationalParameter )
{
    // Define the tolerance of a singularity
    double singularityTolerance = 20.0 * std::numeric_limits< double >::epsilon( );

    // Create USM7 vector and add velocity hodograph elements
    Eigen::Vector7d unifiedStateModelQuaternionsElements;
    unifiedStateModelQuaternionsElements.segment( CHodographQuaternionIndex, 3 ) =
            unifiedStateModelExponentialMapElements.segment( CHodographExponentialMapIndex, 3 );

    // Extract exponential map
    Eigen::Vector3d exponentialMapVector =
            unifiedStateModelExponentialMapElements.segment( e1ExponentialMapIndex, 3 );
    double exponentialMapMagnitude = exponentialMapVector.norm( );

    // Convert exponential map to quaternions
    Eigen::Vector3d epsilonQuaternionVector = Eigen::Vector3d::Zero( );
    if ( std::fabs( exponentialMapMagnitude ) < singularityTolerance )
    {
        epsilonQuaternionVector = exponentialMapVector * ( 0.5 + std::pow( exponentialMapMagnitude, 2 ) / 48.0 );
    }
    {
        epsilonQuaternionVector = exponentialMapVector / exponentialMapMagnitude *
                std::sin( 0.5 * exponentialMapMagnitude );
    }
    double etaQuaternionParameter = std::cos( 0.5 * exponentialMapMagnitude );

    // Add quaternions to USM7 vector
    unifiedStateModelQuaternionsElements.segment( epsilon1QuaternionIndex, 3 ) = epsilonQuaternionVector;
    unifiedStateModelQuaternionsElements( etaQuaternionIndex ) = etaQuaternionParameter;

    // Give back result
    return convertUnifiedStateModelQuaternionsToCartesianElements( unifiedStateModelQuaternionsElements,
                                                                   centralBodyGravitationalParameter );
}

} // close namespace orbital_element_conversions

} // close namespace tudat

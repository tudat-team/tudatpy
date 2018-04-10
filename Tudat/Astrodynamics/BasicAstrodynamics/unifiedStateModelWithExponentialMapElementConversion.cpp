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
#include "Tudat/Astrodynamics/BasicAstrodynamics/unifiedStateModelWithExponentialMapElementConversion.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"

namespace tudat
{

namespace orbital_element_conversions
{

//! Convert Keplerian elements to unified state model elements with exponential map.
Eigen::Vector6d convertKeplerianToUnifiedStateModelWithExponentialMapElements(
        const Eigen::Vector6d& keplerianElements,
        const double centralBodyGravitationalParameter )
{
    using mathematical_constants::PI;

    // Declaring eventual output vector.
    Eigen::Vector6d convertedUnifiedStateModelElements = Eigen::Vector6d::Zero( );

    // Define the tolerance of a singularity
    double singularityTolerance = 1.0e-15; // Based on tolerance chosen in
                                           // orbitalElementConversions.cpp in Tudat Core.

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
    double argumentOfLongitude = keplerianElements( argumentOfPeriapsisIndex ) +
            keplerianElements( trueAnomalyIndex ); // also called u
    double rightAscensionOfLatitude = argumentOfLongitude +
            keplerianElements( longitudeOfAscendingNodeIndex );  // also called lambda
    double XiParameter = ( std::cos( keplerianElements( inclinationIndex ) ) + 1.0 ) *
            ( std::cos( rightAscensionOfLatitude ) + 1 ) - 2.0;

    // Check Xi parameter for numerical errors (magnitude cannot be larger than 1)
    double multiplicativeConstant = 0.0;
    if ( ( std::fabs( std::fabs( XiParameter ) - 2.0 ) < 2.0 * singularityTolerance ) )
    {
        // Remove numerical error and compute multiplicative constant of exponential map
        // multiplicativeConstant = ( XiParameter > 0.0 ) ? 0.5 : Inf;
        if ( XiParameter > 0.0 )
        {
            multiplicativeConstant = 0.5;
            // both numerator and denominator approach zero, and de l'Hopital's rule says the function itself reaches 0.5
        }
        else
        {
            // Define the error message.
            std::stringstream errorMessage;
            errorMessage << "The value of Xi is too close to negative 2."
                         << "This value generates a singularity in the convertion, due to the pure-retrograde nature of the orbit. \n"
                         << "Computed Xi parameter: " << XiParameter << ".\n" << std::endl;

            // Throw exception.
            throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
        }
    }
    else if ( ( std::fabs( XiParameter ) - 2.0 ) <= 0.0 )
    {
        // Compute multiplicative constant of exponential map
        multiplicativeConstant = std::acos( 0.5 * XiParameter ) / std::sqrt( 4.0 - std::pow( XiParameter, 2 ) );
    }
    else
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "The magnitude of the Xi parameter is supposed to be smaller than 2.\n"
                     << "Computed Xi parameter: " << XiParameter << ".\n" << std::endl;

        // Throw exception.
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }

    // Compute the epsilon1 quaternion of the unified state model
    convertedUnifiedStateModelElements( e1ExponentialMapIndex ) =
            multiplicativeConstant * std::sin( keplerianElements( inclinationIndex ) ) *
            ( std::cos( keplerianElements( longitudeOfAscendingNodeIndex ) ) + std::cos( argumentOfLongitude ) );

    // Compute the epsilon2 quaternion of the unified state model
    convertedUnifiedStateModelElements( e2ExponentialMapIndex ) =
            multiplicativeConstant * std::sin( keplerianElements( inclinationIndex ) ) *
            ( std::sin( keplerianElements( longitudeOfAscendingNodeIndex ) ) - std::sin( argumentOfLongitude ) );

    // Compute the epsilon3 quaternion of the unified state model
    convertedUnifiedStateModelElements( e3ExponentialMapIndex ) =
            multiplicativeConstant * std::sin( rightAscensionOfLatitude ) *
            ( std::cos( keplerianElements( inclinationIndex ) ) + 1.0 );

    // Give back result
    return convertedUnifiedStateModelElements;

}

//! Convert unified state model elements with exponential map to Keplerian elements.
Eigen::Vector6d convertUnifiedStateModelWithExponentialMapToKeplerianElements(
        const Eigen::Vector6d& unifiedStateModelElements,
        const double centralBodyGravitationalParameter )
{
    using mathematical_constants::PI;

    // Declaring eventual output vector.
    Eigen::Vector6d convertedKeplerianElements = Eigen::Vector6d::Zero( );

    // Define the tolerance of a singularity
    double singularityTolerance = 20.0 * std::numeric_limits< double >::epsilon( );

    // Compute auxiliary parameters
    double exponentialMapMagnitude = unifiedStateModelElements.segment( 3, 3 ).norm( ); // magnitude of exponential map, also called xi
    double XiParameter = 2.0 * std::cos( exponentialMapMagnitude );

    // Compute right ascension of latitude
    double lambdaParameter = 0.0;
    if ( exponentialMapMagnitude < singularityTolerance )
    {
        lambdaParameter = 2.0 * std::atan( 0.5 * unifiedStateModelElements( e3ExponentialMapIndex ) );
    }
    else
    {
        lambdaParameter = 2.0 * std::atan( unifiedStateModelElements( e3ExponentialMapIndex ) / exponentialMapMagnitude *
                                         std::tan( 0.5 * exponentialMapMagnitude ) );
    }

    // Trigonometric values of right ascension of latitude
    double cosineLambda = std::cos( lambdaParameter );
    double sineLambda = std::sin( lambdaParameter );

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

    // Compute inclination
    if ( ( ( XiParameter - 2.0 ) < singularityTolerance ) || ( ( cosineLambda - 1.0 ) < singularityTolerance ) )
    {
        // Convert to quaternions
        Eigen::Vector3d exponentialMapVector = unifiedStateModelElements.segment( e1ExponentialMapIndex, 3 );
        double exponentialMapMagnitude = exponentialMapVector.norm( ); // also called xi
        Eigen::Vector3d epsilonQuaternionVector = exponentialMapVector / exponentialMapMagnitude *
                std::sin( 0.5 * exponentialMapMagnitude ); // xi is non-zero since Xi is -2.0
        double epsilon1Quaternion = epsilonQuaternionVector( 0 );
        double epsilon2Quaternion = epsilonQuaternionVector( 1 );

        // Compute inclination with quaternions
        convertedKeplerianElements( inclinationIndex ) =
                std::acos( 1.0 - 2.0 * ( epsilon1Quaternion *
                                         epsilon1Quaternion +
                                         epsilon2Quaternion *
                                         epsilon2Quaternion ) );
    }
    else
    {
        convertedKeplerianElements( inclinationIndex ) = std::acos(
                    ( XiParameter + 2.0 ) / ( cosineLambda + 1.0 ) - 1.0 );
        // this acos is always defined correctly because the inclination is always below pi rad.
    }

    // Compute longitude of ascending node
    if ( std::fabs( std::fabs( convertedKeplerianElements( inclinationIndex ) ) - PI ) < singularityTolerance )
        // pure-prograde or pure-retrograde orbit
    {
        convertedKeplerianElements( longitudeOfAscendingNodeIndex ) = 0.0; // by definition
    }
    else
    {
        convertedKeplerianElements( longitudeOfAscendingNodeIndex ) = 0.5 * lambdaParameter +
                std::atan2( unifiedStateModelElements( e2ExponentialMapIndex ),
                            unifiedStateModelElements( e1ExponentialMapIndex ) );

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
            convertedKeplerianElements( longitudeOfAscendingNodeIndex ) =
                    convertedKeplerianElements( longitudeOfAscendingNodeIndex ) + 2.0 * PI;
        }
    }

    // Compute true anomaly and argument of periapsis
    if ( std::fabs( RHodographElement ) < singularityTolerance ) // circular orbit
    {
        convertedKeplerianElements( argumentOfPeriapsisIndex ) = 0.0; // by definition
        convertedKeplerianElements( trueAnomalyIndex ) =
                lambdaParameter - convertedKeplerianElements( longitudeOfAscendingNodeIndex );

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
            convertedKeplerianElements( trueAnomalyIndex ) =
                    convertedKeplerianElements( trueAnomalyIndex ) + 2.0 * PI;
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
            convertedKeplerianElements( trueAnomalyIndex ) =
                    convertedKeplerianElements( trueAnomalyIndex ) + 2.0 * PI;
        }

        convertedKeplerianElements( argumentOfPeriapsisIndex ) =
                lambdaParameter -
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
            convertedKeplerianElements( argumentOfPeriapsisIndex ) =
                    convertedKeplerianElements( argumentOfPeriapsisIndex ) + 2.0 * PI;
        }
    }

    // Give back result
    return convertedKeplerianElements;
}

} // close namespace orbital_element_conversions

} // close namespace tudat

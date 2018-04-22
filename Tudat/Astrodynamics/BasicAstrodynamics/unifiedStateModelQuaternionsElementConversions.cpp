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
 *          and navigation. Master's thesis, Delft University of Technology.
 *      <Second reference>
 *
 */

#include <cmath>
#include <iostream>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/missionGeometry.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unifiedStateModelQuaternionsElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"

namespace tudat
{

namespace orbital_element_conversions
{

//! Convert Keplerian elements to unified state model elements with quaternions.
Eigen::Vector7d convertKeplerianToUnifiedStateModelQuaternionsElements(
        const Eigen::Vector6d& keplerianElements,
        const double centralBodyGravitationalParameter )
{
    using mathematical_constants::PI;

    // Declaring eventual output vector.
    Eigen::Vector7d convertedUnifiedStateModelElements = Eigen::Vector7d::Zero( );

    // Define the tolerance of a singularity
    const double singularityTolerance = 20.0 * std::numeric_limits< double >::epsilon( );

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
        convertedUnifiedStateModelElements( CHodographQuaternionIndex ) =
                std::sqrt( centralBodyGravitationalParameter / keplerianElements( semiLatusRectumIndex ) );
    }
    else
    {
        convertedUnifiedStateModelElements( CHodographQuaternionIndex ) =
                std::sqrt( centralBodyGravitationalParameter / ( keplerianElements( semiMajorAxisIndex )
                                                                 * ( 1 - keplerianElements( eccentricityIndex ) *
                                                                     keplerianElements( eccentricityIndex ) ) ) );
    }

    // Calculate the additional R hodograph parameter
    double RHodographElement = keplerianElements( eccentricityIndex ) *
            convertedUnifiedStateModelElements( CHodographQuaternionIndex );

    // Compute the Rf1 hodograph element of the unified state model
    convertedUnifiedStateModelElements( Rf1HodographQuaternionIndex ) =
            - RHodographElement * std::sin( keplerianElements( longitudeOfAscendingNodeIndex )
                                            + keplerianElements( argumentOfPeriapsisIndex ) );

    // Compute the Rf2 hodograph element of the unified state model
    convertedUnifiedStateModelElements( Rf2HodographQuaternionIndex ) =
            RHodographElement * std::cos( keplerianElements( longitudeOfAscendingNodeIndex )
                                          + keplerianElements( argumentOfPeriapsisIndex ) );

    // Calculate the additional argument of longitude u
    double argumentOfLatitude = keplerianElements( argumentOfPeriapsisIndex ) +
            keplerianElements( trueAnomalyIndex );

    // Compute the epsilon1 quaternion of the unified state model
    convertedUnifiedStateModelElements( epsilon1QuaternionIndex ) =
            std::sin( 0.5 * keplerianElements( inclinationIndex ) ) *
            std::cos( 0.5 * ( keplerianElements( longitudeOfAscendingNodeIndex ) - argumentOfLatitude ) );

    // Compute the epsilon2 quaternion of the unified state model
    convertedUnifiedStateModelElements( epsilon2QuaternionIndex ) =
            std::sin( 0.5 * keplerianElements( inclinationIndex ) ) *
            std::sin( 0.5 * ( keplerianElements( longitudeOfAscendingNodeIndex ) - argumentOfLatitude ) );

    // Compute the epsilon3 quaternion of the unified state model
    convertedUnifiedStateModelElements( epsilon3QuaternionIndex ) =
            std::cos( 0.5 * keplerianElements( inclinationIndex ) ) *
            std::sin( 0.5 * ( keplerianElements( longitudeOfAscendingNodeIndex ) + argumentOfLatitude ) );

    // Compute the eta quaternion of the unified state model
    convertedUnifiedStateModelElements( etaQuaternionIndex ) =
            std::cos( 0.5 * keplerianElements( inclinationIndex ) ) *
            std::cos( 0.5 * ( keplerianElements( longitudeOfAscendingNodeIndex ) + argumentOfLatitude ) );

    // Give back result
    return convertedUnifiedStateModelElements;
}

//! Convert unified state model elements with quaternions to Keplerian elements.
Eigen::Vector6d convertUnifiedStateModelQuaternionsToKeplerianElements(
        const Eigen::Vector7d& unifiedStateModelElements,
        const double centralBodyGravitationalParameter,
        const bool forceQuaternionNormalization )
{
    using mathematical_constants::PI;

    // Declaring eventual output vector.
    Eigen::Vector6d convertedKeplerianElements = Eigen::Vector6d::Zero( );

    // Define the tolerance of a singularity
    const double singularityTolerance = 20.0 * std::numeric_limits< double >::epsilon( );

    // Check whether the unified state model elements are within expected limits
    // If inclination is zero and the right ascension of ascending node is non-zero
    Eigen::Vector4d quaternionsVector = unifiedStateModelElements.segment( epsilon1QuaternionIndex, 4 );
    const double normOfQuaternionElements = quaternionsVector.norm( );
    if ( std::fabs( normOfQuaternionElements - 1.0 ) > singularityTolerance )
    {
        if ( forceQuaternionNormalization )
        {
            quaternionsVector /= normOfQuaternionElements;
        }
        else
        {
            // Define the error message.
            std::stringstream errorMessage;
            errorMessage << "The norm of the quaternion should be equal to one.\n"
                         << "Norm of the specified quaternion is: " << normOfQuaternionElements - 1.0 << " + 1." << std::endl;

            // Throw exception.
            throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
        }
    }
    // Else, nothing wrong and continue

    // Extract quaternion elements
    double epsilon1QuaternionParameter = quaternionsVector( 0 );
    double epsilon2QuaternionParameter = quaternionsVector( 1 );
    double epsilon3QuaternionParameter = quaternionsVector( 2 );
    double etaQuaternionParameter = quaternionsVector( 3 );

    // Check whether the orbit is pure-retrograde
    if ( ( std::fabs( epsilon3QuaternionParameter ) < singularityTolerance ) &&
         ( std::fabs( etaQuaternionParameter ) < singularityTolerance ) )
        // pure-retrograde orbit -> inclination = PI
    {
        // Define the error message
        std::stringstream errorMessage;
        errorMessage << "Pure-retrograde orbit (inclination = PI).\n"
                     << "Unified state model elements cannot be transformed to Kepler elements." << std::endl;

        // Throw exception
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }

    // Compute auxiliary parameters cosineLambda, sineLambda and lambda
    double denominator = epsilon3QuaternionParameter * epsilon3QuaternionParameter +
            etaQuaternionParameter * etaQuaternionParameter;
    double cosineLambda = ( etaQuaternionParameter * etaQuaternionParameter -
                            epsilon3QuaternionParameter * epsilon3QuaternionParameter ) / denominator;
    double sineLambda = ( 2.0 * epsilon3QuaternionParameter *
                          etaQuaternionParameter ) / denominator;
    double rightAscensionOfLatitude = std::atan2( sineLambda, cosineLambda );

    // Compute auxiliary parameters auxiliaryParameter1 and auxiliaryParameter2
    double auxiliaryParameter1 = unifiedStateModelElements( Rf1HodographQuaternionIndex ) * cosineLambda +
            unifiedStateModelElements( Rf2HodographQuaternionIndex ) * sineLambda;
    double auxiliaryParameter2 = unifiedStateModelElements( CHodographQuaternionIndex ) -
            unifiedStateModelElements( Rf1HodographQuaternionIndex ) * sineLambda +
            unifiedStateModelElements( Rf2HodographQuaternionIndex ) * cosineLambda;

    // Compute auxiliary R hodograph parameter
    double RHodographElement = std::sqrt( unifiedStateModelElements( Rf1HodographQuaternionIndex ) *
                                          unifiedStateModelElements( Rf1HodographQuaternionIndex ) +
                                          unifiedStateModelElements( Rf2HodographQuaternionIndex ) *
                                          unifiedStateModelElements( Rf2HodographQuaternionIndex ) );

    // Compute eccentricity
    convertedKeplerianElements( eccentricityIndex ) =
            RHodographElement / unifiedStateModelElements( CHodographQuaternionIndex );

    // Compute semi-major axis or, in case of a parabolic orbit, the semi-latus rectum.
    if ( std::fabs( convertedKeplerianElements( eccentricityIndex ) - 1.0 ) < singularityTolerance )
        // parabolic orbit -> semi-major axis is not defined. Use semi-latus rectum instead.
    {
        convertedKeplerianElements( semiLatusRectumIndex ) = centralBodyGravitationalParameter /
                ( unifiedStateModelElements( CHodographQuaternionIndex ) * unifiedStateModelElements( CHodographQuaternionIndex ) );
    }
    else
    {
        convertedKeplerianElements( semiMajorAxisIndex ) =
                centralBodyGravitationalParameter /
                ( std::pow( unifiedStateModelElements( CHodographQuaternionIndex ), 2 ) *
                  ( 1 - std::pow( convertedKeplerianElements( eccentricityIndex ), 2 ) ) );
    }

    // Compute inclination
    convertedKeplerianElements( inclinationIndex ) =
            std::acos( 1.0 - 2.0 * ( epsilon1QuaternionParameter * epsilon1QuaternionParameter +
                                     epsilon2QuaternionParameter * epsilon2QuaternionParameter ) );
    // This acos is always defined correctly because the inclination is always below pi rad.

    // Find sine and cosine of longitude of ascending node separately
    double sineOmega = epsilon1QuaternionParameter * epsilon3QuaternionParameter +
            epsilon2QuaternionParameter * etaQuaternionParameter;
    double cosineOmega = epsilon1QuaternionParameter * etaQuaternionParameter -
            epsilon2QuaternionParameter * epsilon3QuaternionParameter;
    denominator = std::sqrt( cosineOmega * cosineOmega +
                             sineOmega * sineOmega ); // overwrite old denominator

    // Compute longitude of ascending node
    if ( std::fabs( convertedKeplerianElements( inclinationIndex ) - PI ) < singularityTolerance )
        // pure-retrograde orbit -> inclination = PI
    {
        // Define the error message
        std::stringstream errorMessage;
        errorMessage << "Pure-retrograde orbit (inclination = PI).\n"
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
                            ( ( auxiliaryParameter2 - unifiedStateModelElements( CHodographQuaternionIndex ) )
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

////! Convert Cartesian elements to unified state model elements with quaternions.
Eigen::Vector7d convertCartesianToUnifiedStateModelQuaternionsElements(
        const Eigen::Vector6d& cartesianElements,
        const double centralBodyGravitationalParameter )
{
    using mathematical_constants::PI;

    // Declaring eventual output vector.
    Eigen::Vector7d convertedUnifiedStateModelElements = Eigen::Vector7d::Zero( );

    // Define the tolerance of a singularity
    const double singularityTolerance = 20.0 * std::numeric_limits< double >::epsilon( );

    // Find Cartesian position and velocity vectors and magnitudes
    Eigen::Vector3d positionVector = cartesianElements.segment( xCartesianPositionIndex, 3 );
    double positionMagnitude = positionVector.norm( );
    Eigen::Vector3d velocityVector = cartesianElements.segment( xCartesianVelocityIndex, 3 );

    // Determine specific angular momentum vector and magnitude
    Eigen::Vector3d angularMomentumVector = positionVector.cross( velocityVector );
    double angularMomentumMagnitude = angularMomentumVector.norm( );

    // Check whether the orbit is pure-retrograde
    double angularMomentumMagnitudeMinusZComponent = angularMomentumMagnitude + angularMomentumVector( 2 );
    if ( std::fabs( angularMomentumMagnitudeMinusZComponent ) < singularityTolerance )
        // pure-retrograde orbit -> inclination = PI
    {
        // Define the error message
        std::stringstream errorMessage;
        errorMessage << "Pure-retrograde orbit (inclination = PI).\n"
                     << "Unified state model elements cannot be transformed to Kepler elements." << std::endl;

        // Throw exception
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }

    // Find C hodograph element of the unified state model
    convertedUnifiedStateModelElements( CHodographQuaternionIndex ) = centralBodyGravitationalParameter /
            angularMomentumMagnitude;

    // Find direction cosine matrix with position and angular momentum vectors
    Eigen::Matrix3d directionCosineMatrix = Eigen::Matrix3d::Zero( );
    directionCosineMatrix.block( 0, 0, 1, 3 ) = positionVector.normalized( ).transpose( );
    directionCosineMatrix.block( 1, 0, 1, 3 ) = ( angularMomentumVector.normalized( ).cross(
                                                      positionVector.normalized( ) ) ).transpose( );
    directionCosineMatrix.block( 2, 0, 1, 3 ) = angularMomentumVector.normalized( ).transpose( );

    // Compute square of quaternions
    double traceDirectionCosineMatrix = directionCosineMatrix.trace( );
    Eigen::Vector4d quaternionSquaredVector;
    for ( unsigned int i = 0; i < 3; i++ )
    {
        quaternionSquaredVector( i ) = ( 1.0 - traceDirectionCosineMatrix + 2.0 *
                                         directionCosineMatrix( i, i ) ) / 4.0;
    }
    quaternionSquaredVector( 3 ) = ( 1.0 + traceDirectionCosineMatrix ) / 4.0;

    // Based on the maximum value, find the quaternion elements
    Eigen::Vector4d::Index indexLargestQuaternionElement;
    double valueLargestQuaternionElement = quaternionSquaredVector.maxCoeff( &indexLargestQuaternionElement );
    switch ( indexLargestQuaternionElement )
    {
    case 0:
    {
        // Find value of largest quaternion parameter
        convertedUnifiedStateModelElements( epsilon1QuaternionIndex ) = std::sqrt( valueLargestQuaternionElement );

        // Find other values
        Eigen::Vector3d auxiliaryVector = Eigen::Vector3d::Zero( );
        auxiliaryVector( 0 ) = directionCosineMatrix( 1, 0 ) + directionCosineMatrix( 0, 1 );
        auxiliaryVector( 1 ) = directionCosineMatrix( 2, 0 ) + directionCosineMatrix( 0, 2 );
        auxiliaryVector( 2 ) = directionCosineMatrix( 1, 2 ) - directionCosineMatrix( 2, 1 );
        auxiliaryVector /= 4 * convertedUnifiedStateModelElements( epsilon1QuaternionIndex );

        // Distribute to state vector
        convertedUnifiedStateModelElements( epsilon2QuaternionIndex ) = auxiliaryVector( 0 );
        convertedUnifiedStateModelElements( epsilon3QuaternionIndex ) = auxiliaryVector( 1 );
        convertedUnifiedStateModelElements( etaQuaternionIndex ) = auxiliaryVector( 2 );
        break;
    }
    case 1:
    {
        // Find value of largest quaternion parameter
        convertedUnifiedStateModelElements( epsilon2QuaternionIndex ) = std::sqrt( valueLargestQuaternionElement );

        // Find other values
        Eigen::Vector3d auxiliaryVector = Eigen::Vector3d::Zero( );
        auxiliaryVector( 0 ) = directionCosineMatrix( 0, 1 ) + directionCosineMatrix( 1, 0 );
        auxiliaryVector( 1 ) = directionCosineMatrix( 2, 1 ) + directionCosineMatrix( 1, 2 );
        auxiliaryVector( 2 ) = directionCosineMatrix( 2, 0 ) - directionCosineMatrix( 0, 2 );
        auxiliaryVector /= 4 * convertedUnifiedStateModelElements( epsilon2QuaternionIndex );

        // Distribute to state vector
        convertedUnifiedStateModelElements( epsilon1QuaternionIndex ) = auxiliaryVector( 0 );
        convertedUnifiedStateModelElements( epsilon3QuaternionIndex ) = auxiliaryVector( 1 );
        convertedUnifiedStateModelElements( etaQuaternionIndex ) = auxiliaryVector( 2 );
        break;
    }
    case 2:
    {
        // Find value of largest quaternion parameter
        convertedUnifiedStateModelElements( epsilon3QuaternionIndex ) = std::sqrt( valueLargestQuaternionElement );

        // Find other values
        Eigen::Vector3d auxiliaryVector = Eigen::Vector3d::Zero( );
        auxiliaryVector( 0 ) = directionCosineMatrix( 0, 2 ) + directionCosineMatrix( 2, 0 );
        auxiliaryVector( 1 ) = directionCosineMatrix( 1, 2 ) + directionCosineMatrix( 2, 1 );
        auxiliaryVector( 2 ) = directionCosineMatrix( 0, 1 ) - directionCosineMatrix( 1, 0 );
        auxiliaryVector /= 4 * convertedUnifiedStateModelElements( epsilon3QuaternionIndex );

        // Distribute to state vector
        convertedUnifiedStateModelElements( epsilon1QuaternionIndex ) = auxiliaryVector( 0 );
        convertedUnifiedStateModelElements( epsilon2QuaternionIndex ) = auxiliaryVector( 1 );
        convertedUnifiedStateModelElements( etaQuaternionIndex ) = auxiliaryVector( 2 );
        break;
    }
    case 3:
    {
        // Find value of largest quaternion parameter
        convertedUnifiedStateModelElements( etaQuaternionIndex ) = std::sqrt( valueLargestQuaternionElement );

        // Find other values
        Eigen::Vector3d auxiliaryVector = Eigen::Vector3d::Zero( );
        auxiliaryVector( 0 ) = directionCosineMatrix( 1, 2 ) - directionCosineMatrix( 2, 1 );
        auxiliaryVector( 1 ) = directionCosineMatrix( 2, 0 ) - directionCosineMatrix( 0, 2 );
        auxiliaryVector( 2 ) = directionCosineMatrix( 0, 1 ) - directionCosineMatrix( 1, 0 );
        auxiliaryVector /= 4 * convertedUnifiedStateModelElements( etaQuaternionIndex );

        // Distribute to state vector
        convertedUnifiedStateModelElements( epsilon1QuaternionIndex ) = auxiliaryVector( 0 );
        convertedUnifiedStateModelElements( epsilon2QuaternionIndex ) = auxiliaryVector( 1 );
        convertedUnifiedStateModelElements( epsilon3QuaternionIndex ) = auxiliaryVector( 2 );
        break;
    }
    default:
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "Could not find the maximum value of the quaternion.\n"
                     << "Specified squared quaternion: " << quaternionSquaredVector.transpose( ) << std::endl;

        // Throw exception.
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }
    }

    // Recompute sine and cosine of right ascension of latitude (lambda), for better numerical accuracy
    double denominator = convertedUnifiedStateModelElements( epsilon3QuaternionIndex ) *
            convertedUnifiedStateModelElements( epsilon3QuaternionIndex ) +
            convertedUnifiedStateModelElements( etaQuaternionIndex ) * convertedUnifiedStateModelElements( etaQuaternionIndex );
    double cosineLambda = ( convertedUnifiedStateModelElements( etaQuaternionIndex ) *
                            convertedUnifiedStateModelElements( etaQuaternionIndex ) -
                            convertedUnifiedStateModelElements( epsilon3QuaternionIndex ) *
                            convertedUnifiedStateModelElements( epsilon3QuaternionIndex ) ) / denominator;
    double sineLambda = ( 2.0 * convertedUnifiedStateModelElements( epsilon3QuaternionIndex ) *
                          convertedUnifiedStateModelElements( etaQuaternionIndex ) ) / denominator;

    // Compute auxiliary parameters
    double radialVelocity = positionVector.dot( velocityVector ) / positionMagnitude;
    Eigen::Vector3d auxiliaryParameter3 = radialVelocity / positionMagnitude * positionVector;
    double auxiliaryParameter2 = ( velocityVector - auxiliaryParameter3 ).norm( );
    double auxiliaryParameter1 = std::signbit( radialVelocity ) ?
                - auxiliaryParameter3.norm( ) : auxiliaryParameter3.norm( ); // take norm now that vector value has been used
    // The sign of first velocity component depends on true anomaly (positive if < PI), and true anomaly can be related to the radial
    // velocity

    // Compute Rf1 and Rf2 hodograph elements
    convertedUnifiedStateModelElements( Rf1HodographQuaternionIndex ) = auxiliaryParameter1 * cosineLambda -
            ( auxiliaryParameter2 - convertedUnifiedStateModelElements( CHodographQuaternionIndex ) ) * sineLambda;
    convertedUnifiedStateModelElements( Rf2HodographQuaternionIndex ) = auxiliaryParameter1 * sineLambda +
            ( auxiliaryParameter2 - convertedUnifiedStateModelElements( CHodographQuaternionIndex ) ) * cosineLambda;

    // Give back result
    return convertedUnifiedStateModelElements;
}

//! Convert unified state model elements with quaternions to Cartesian elements.
Eigen::Vector6d convertUnifiedStateModelQuaternionsToCartesianElements(
        const Eigen::Vector7d& unifiedStateModelElements,
        const double centralBodyGravitationalParameter,
        const bool forceQuaternionNormalization )
{
    using mathematical_constants::PI;

    // Declaring eventual output vector.
    Eigen::Vector6d convertedCartesianElements = Eigen::Vector6d::Zero( );

    // Define the tolerance of a singularity
    const double singularityTolerance = 20.0 * std::numeric_limits< double >::epsilon( );

    // Extract quaternion elements
    Eigen::Vector4d quaternionsVector = unifiedStateModelElements.segment( epsilon1QuaternionIndex, 4 );
    const double normOfQuaternionElements = quaternionsVector.norm( );
    if ( std::fabs( normOfQuaternionElements - 1.0 ) > singularityTolerance )
    {
        if ( forceQuaternionNormalization )
        {
            quaternionsVector /= normOfQuaternionElements;
        }
        else
        {
            // Define the error message.
            std::stringstream errorMessage;
            errorMessage << "The norm of the quaternion should be equal to one.\n"
                         << "Norm of the specified quaternion is: " << normOfQuaternionElements - 1.0 << " + 1." << std::endl;

            // Throw exception.
            throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
        }
    }
    // Else, nothing wrong and continue

    // Extract quaternion elements
    double epsilon3QuaternionParameter = quaternionsVector( 2 );
    double etaQuaternionParameter = quaternionsVector( 3 );

    // Declare auxiliary parameters before using them in the if statement
    double cosineLambda = 0.0;
    double sineLambda = 0.0;
    double rightAscensionOfLatitude = 0.0;

    // Compute auxiliary parameters cosineLambda, sineLambda and lambda
    if ( ( std::fabs( epsilon3QuaternionParameter ) < singularityTolerance ) &&
         ( std::fabs( etaQuaternionParameter ) < singularityTolerance ) )
        // pure-retrograde orbit -> inclination = PI
    {
        // Define the error message
        std::stringstream errorMessage;
        errorMessage << "Pure-retrograde orbit (inclination = PI).\n"
                     << "Unified state model elements cannot be transformed to Kepler elements." << std::endl;

        // Throw exception
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }
    else
    {
        double denominator = epsilon3QuaternionParameter * epsilon3QuaternionParameter +
                etaQuaternionParameter * etaQuaternionParameter;
        cosineLambda = ( etaQuaternionParameter * etaQuaternionParameter -
                         epsilon3QuaternionParameter * epsilon3QuaternionParameter ) / denominator;
        sineLambda = ( 2.0 * epsilon3QuaternionParameter *
                       etaQuaternionParameter ) / denominator;
        rightAscensionOfLatitude = std::atan2( sineLambda, cosineLambda );
    }

    // Compute auxiliary parameters auxiliaryParameter1, auxiliaryParameter2 and auxiliaryVector1
    double auxiliaryParameter1 = unifiedStateModelElements( Rf1HodographQuaternionIndex ) * cosineLambda +
            unifiedStateModelElements( Rf2HodographQuaternionIndex ) * sineLambda;
    double auxiliaryParameter2 = unifiedStateModelElements( CHodographQuaternionIndex ) -
            unifiedStateModelElements( Rf1HodographQuaternionIndex ) * sineLambda +
            unifiedStateModelElements( Rf2HodographQuaternionIndex ) * cosineLambda;
    Eigen::Vector2d auxiliaryVector1;
    auxiliaryVector1( 0 ) = auxiliaryParameter1;
    auxiliaryVector1( 1 ) = auxiliaryParameter2;

    // Find direction cosine matrix in terms of quaternions
    Eigen::Matrix3d inverseDirectionCosineMatrix =
            linear_algebra::computeDirectionCosineMatrixFromQuaternions( quaternionsVector, true );

    // Get Cartesian position vector
    convertedCartesianElements.segment( xCartesianPositionIndex, 3 ) =
            centralBodyGravitationalParameter / unifiedStateModelElements( CHodographQuaternionIndex ) /
            auxiliaryParameter2 * inverseDirectionCosineMatrix.block( 0, 0, 3, 1 ); // take first column of matrix

    // Get Cartesian velocity vector
    convertedCartesianElements.segment( xCartesianVelocityIndex, 3 ) =
            inverseDirectionCosineMatrix.block( 0, 0, 3, 2 ) * auxiliaryVector1;

    // Give back result
    return convertedCartesianElements;
}

} // close namespace orbital_element_conversions

} // close namespace tudat

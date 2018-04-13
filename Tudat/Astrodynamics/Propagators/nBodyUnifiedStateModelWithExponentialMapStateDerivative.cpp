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
 */

#include "Tudat/Astrodynamics/Propagators/nBodyUnifiedStateModelWithExponentialMapStateDerivative.h"

namespace tudat
{

namespace propagators
{

//! Function to evaluate the state derivative for the unified state model with exponential map
Eigen::Vector6d computeStateDerivativeForUnifiedStateModelWithExponentialMap(
        const Eigen::Vector6d& currentUnifiedStateModelElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double sineLambdaParameter,
        const double cosineLambdaParameter,
        const double gammaParameter,
        const Eigen::Vector3d& rotationalVelocityVector,
        const Eigen::Vector3d& pParameterVector )
{
//    // REMOVE vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv REMOVE
//    std::cout << "USM state: " << std::endl << currentUnifiedStateModelElements << std::endl;
//    std::cout << "Acceleration: " << std::endl << accelerationsInRswFrame << std::endl;
//    std::cout << "Sine lambda: " << sineLambdaParameter << std::endl;
//    std::cout << "Cosine lambda: " << cosineLambdaParameter << std::endl;
//    std::cout << "Gamma: " << gammaParameter << std::endl;
//    std::cout << "Rotational velocity: " << std::endl << rotationalVelocityVector << std::endl;
//    std::cout << "P parameter: " << std::endl << pParameterVector << std::endl;
//    // REMOVE ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ REMOVE

    // Define the tolerance of a singularity
    double singularityTolerance = 20.0 * std::numeric_limits< double >::epsilon( );
//    // REMOVE vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv REMOVE
//    std::cout << "Tolerance: " << std::endl << singularityTolerance << std::endl;
//    // REMOVE ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ REMOVE

    // Compute matrix for dynamic equation
    Eigen::Matrix3d hodographMatrix = Eigen::Matrix3d::Zero( );
    hodographMatrix( 0, 1 ) = - pParameterVector( 0 );
    hodographMatrix( 1, 0 ) =   cosineLambdaParameter;
    hodographMatrix( 1, 1 ) = - ( 1 + pParameterVector( 0 ) ) * sineLambdaParameter;
    hodographMatrix( 1, 2 ) = - gammaParameter * pParameterVector( 1 );
    hodographMatrix( 2, 0 ) =   sineLambdaParameter;
    hodographMatrix( 2, 1 ) =   ( 1 + pParameterVector( 0 ) ) * cosineLambdaParameter;
    hodographMatrix( 2, 2 ) =   gammaParameter * pParameterVector( 2 );
//    // REMOVE vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv REMOVE
//    std::cout << "Hodograph matrix: " << std::endl << hodographMatrix << std::endl;
//    // REMOVE ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ REMOVE

    // Compute kinematic equation, i.e., derivative of exponential map
    Eigen::Vector3d exponentialMapVector = currentUnifiedStateModelElements.segment( 3, 3 );
    double exponentialMapMagnitude = exponentialMapVector.norm( );
    Eigen::Vector3d exponentialMapDerivative = Eigen::Vector3d::Zero( );
    if ( std::fabs( exponentialMapMagnitude ) < singularityTolerance )
    {
        double exponentialMapMagnitudeSquared = std::pow( exponentialMapMagnitude, 2 );
        exponentialMapDerivative = 0.5 * ( ( ( 12.0 - exponentialMapMagnitudeSquared ) / 6.0 ) *
                                           rotationalVelocityVector - rotationalVelocityVector.cross( exponentialMapVector ) -
                                           rotationalVelocityVector.dot( exponentialMapVector ) *
                                           ( ( 60.0 + exponentialMapMagnitudeSquared ) / 360.0 ) * exponentialMapVector );
    }
    else
    {
        double cotangentHalfExponentialMapMagnitude = std::cos( 0.5 * exponentialMapMagnitude ) /
                std::sin( 0.5 * exponentialMapMagnitude );
        Eigen::Matrix3d skewExponentialMapVector = linear_algebra::getCrossProductMatrix( exponentialMapVector );
        exponentialMapDerivative = ( Eigen::Matrix3d::Identity( ) + 0.5 * skewExponentialMapVector +
                                     ( 1 - 0.5 * exponentialMapMagnitude * cotangentHalfExponentialMapMagnitude ) /
                                     ( exponentialMapMagnitude * exponentialMapMagnitude ) * skewExponentialMapVector *
                                     skewExponentialMapVector ) * rotationalVelocityVector;
//        Eigen::Vector3d exponentialMapCrossRotationalVelocityVector = exponentialMapVector.cross( rotationalVelocityVector );
//        exponentialMapDerivative = rotationalVelocityVector + 0.5 * exponentialMapCrossRotationalVelocityVector +
//                ( 1 - 0.5 * exponentialMapMagnitude * cotangentHalfExponentialMapMagnitude ) /
//                std::pow( exponentialMapMagnitude, 2 ) *
//                exponentialMapVector.cross( exponentialMapCrossRotationalVelocityVector );
    }
//    // REMOVE vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv REMOVE
//    std::cout << "Exponential map vector: " << std::endl << exponentialMapVector << std::endl;
//    std::cout << "Exponential map magnitude: " << std::endl << exponentialMapMagnitude << std::endl;
//    std::cout << "Exponential map derivative: " << std::endl << exponentialMapDerivative << std::endl;
//    // REMOVE ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ REMOVE

    // Evaluate USMEM equations.
    Eigen::Vector6d stateDerivative;
    stateDerivative.segment( 0, 3 ) = hodographMatrix * accelerationsInRswFrame;
    stateDerivative.segment( 3, 3 ) = exponentialMapDerivative;
//    // REMOVE vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv REMOVE
//    std::cout << "State derivative: " << std::endl << stateDerivative << std::endl;
//    // REMOVE ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ REMOVE

    return stateDerivative;
}

//! Function to evaluate the state derivative for the unified state model with exponential map
Eigen::Vector6d computeStateDerivativeForUnifiedStateModelWithExponentialMap(
        const Eigen::Vector6d& currentUnifiedStateModelElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double centralBodyGravitationalParameter )
{
    using namespace orbital_element_conversions;

    // Define the tolerance of a singularity
    double singularityTolerance = 20.0 * std::numeric_limits< double >::epsilon( );

    // Retrieve USM elements
    double CHodograph = currentUnifiedStateModelElements( CHodographExponentialMapIndex );
    double Rf1Hodograph = currentUnifiedStateModelElements( Rf1HodographExponentialMapIndex );
    double Rf2Hodograph = currentUnifiedStateModelElements( Rf2HodographExponentialMapIndex );
    Eigen::Vector3d exponentialMapVector = currentUnifiedStateModelElements.segment( e1ExponentialMapIndex, 3 );
    double exponentialMapMagnitude = exponentialMapVector.norm( ); // also called xi

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
    double epsilon1Quaternion = epsilonQuaternionVector( 0 );
    double epsilon2Quaternion = epsilonQuaternionVector( 1 );
    double epsilon3Quaternion = epsilonQuaternionVector( 2 );
    double etaQuaternion = std::cos( 0.5 * exponentialMapMagnitude );

    // Compute supporting parameters
    double quaterionParameter = std::pow( epsilon3Quaternion, 2) + std::pow( etaQuaternion, 2 );
    double sineLambdaParameter = ( 2 * epsilon3Quaternion * etaQuaternion ) / quaterionParameter;
    double cosineLambdaParameter =  ( std::pow( etaQuaternion, 2 ) - std::pow( epsilon3Quaternion, 2 ) ) /
            quaterionParameter;
    double gammaParameter = ( epsilon1Quaternion * epsilon3Quaternion -
                              epsilon2Quaternion * etaQuaternion ) / quaterionParameter;

    double velocityHodographParameter = CHodograph - Rf1Hodograph * sineLambdaParameter +
            Rf2Hodograph * cosineLambdaParameter;
    Eigen::Vector3d rotationalVelocityVector = Eigen::Vector3d::Zero( );
    rotationalVelocityVector( 0 ) = accelerationsInRswFrame( 2 ) / velocityHodographParameter;
    rotationalVelocityVector( 2 ) = std::pow( velocityHodographParameter, 2 ) * CHodograph /
            centralBodyGravitationalParameter;

    Eigen::Vector3d pParameterVector = Eigen::Vector3d::Zero( );
    pParameterVector( 0 ) = CHodograph;
    pParameterVector( 1 ) = Rf2Hodograph;
    pParameterVector( 2 ) = Rf1Hodograph;
    pParameterVector = pParameterVector / velocityHodographParameter;

    // Evaluate USMEM equations
    return computeStateDerivativeForUnifiedStateModelWithExponentialMap(
                currentUnifiedStateModelElements, accelerationsInRswFrame, sineLambdaParameter,
                cosineLambdaParameter, gammaParameter, rotationalVelocityVector, pParameterVector );
}

//! Function to evaluate the state derivative for the unified state model with exponential map
Eigen::Vector6d computeStateDerivativeForUnifiedStateModelWithExponentialMap(
        const Eigen::Vector6d& currentUnifiedStateModelElements,
        const Eigen::Vector6d& currentCartesianState,
        const Eigen::Vector3d& accelerationsInInertialFrame,
        const double centralBodyGravitationalParameter )
{
    return computeStateDerivativeForUnifiedStateModelWithExponentialMap(
                currentUnifiedStateModelElements,
                reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrx(
                    currentCartesianState ) * accelerationsInInertialFrame, centralBodyGravitationalParameter );
}


} // namespace propagators

} // namespace tudat

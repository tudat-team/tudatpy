/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/Propagators/nBodyUnifiedStateModelWithQuaternionsStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/rotationalMotionStateDerivative.h"

namespace tudat
{

namespace propagators
{

//! Function to evaluate the state derivative for the unified state model with quaternions
Eigen::Vector7d computeStateDerivativeForUnifiedStateModelWithQuaternions(
        const Eigen::Vector7d& currentUnifiedStateModelElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double sineLambdaParameter,
        const double cosineLambdaParameter,
        const double gammaParameter,
        const Eigen::Vector3d& rotationalVelocity,
        const Eigen::Vector3d& pParameter )
{
    // REMOVE vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv REMOVE
    std::cout << "USM state: " << std::endl << currentUnifiedStateModelElements << std::endl;
    std::cout << "Norm before propagation: " << currentUnifiedStateModelElements.segment( 3, 4 ).norm( ) << std::endl;
    std::cout << "Acceleration: " << std::endl << accelerationsInRswFrame << std::endl;
    std::cout << "Sine lambda: " << sineLambdaParameter << std::endl;
    std::cout << "Cosine lambda: " << cosineLambdaParameter << std::endl;
    std::cout << "Gamma: " << gammaParameter << std::endl;
    std::cout << "Rotational velocity: " << std::endl << rotationalVelocity << std::endl;
    std::cout << "P parameter: " << std::endl << pParameter << std::endl;
    // REMOVE ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ REMOVE

    // Compute supporting parameters
    Eigen::Matrix3d hodographMatrix = Eigen::Matrix3d::Zero( );
    hodographMatrix( 0, 1 ) = - pParameter( 0 );
    hodographMatrix( 1, 0 ) = cosineLambdaParameter;
    hodographMatrix( 1, 1 ) = - ( 1 + pParameter( 0 ) ) * sineLambdaParameter;
    hodographMatrix( 1, 2 ) = - gammaParameter * pParameter( 1 );
    hodographMatrix( 2, 0 ) = sineLambdaParameter;
    hodographMatrix( 2, 1 ) = ( 1 + pParameter( 0 ) ) * cosineLambdaParameter;
    hodographMatrix( 2, 1 ) = gammaParameter * pParameter( 2 );

    Eigen::MatrixXd quaternionMatrix = getQuaterionToQuaternionRateMatrix( rotationalVelocity );
//    Eigen::MatrixXd::Zero( 4, 4 );
//    quaternionMatrix( 0, 1 ) = rotationalVelocity( 2 );
//    quaternionMatrix( 0, 2 ) = - rotationalVelocity( 1 );
//    quaternionMatrix( 0, 3 ) = rotationalVelocity( 0 );
//    quaternionMatrix( 1, 0 ) = - rotationalVelocity( 2 );
//    quaternionMatrix( 1, 2 ) = rotationalVelocity( 0 );
//    quaternionMatrix( 1, 3 ) = rotationalVelocity( 1 );
//    quaternionMatrix( 2, 0 ) = rotationalVelocity( 1 );
//    quaternionMatrix( 2, 1 ) = - rotationalVelocity( 0 );
//    quaternionMatrix( 2, 3 ) = rotationalVelocity( 2 );
//    quaternionMatrix( 3, 0 ) = - rotationalVelocity( 0 );
//    quaternionMatrix( 3, 1 ) = - rotationalVelocity( 1 );
//    quaternionMatrix( 3, 2 ) = - rotationalVelocity( 2 );

    // Evaluate USM7 equations.
    Eigen::Vector7d stateDerivative;
    stateDerivative.segment( 0, 3 ) = hodographMatrix * accelerationsInRswFrame;
    stateDerivative.segment( 3, 4 ) = 0.5 * quaternionMatrix * currentUnifiedStateModelElements.segment( 3, 4 );

    // Give output
    return stateDerivative;
}

//! Function to evaluate the state derivative for the unified state model with quaternions
Eigen::Vector7d computeStateDerivativeForUnifiedStateModelWithQuaternions(
        const Eigen::Vector7d& currentUnifiedStateModelElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double centralBodyGravitationalParameter )
{
    using namespace orbital_element_conversions;

    // Retrieve USM elements
    double CHodograph = currentUnifiedStateModelElements( CHodographQuaternionIndex );
    double Rf1Hodograph = currentUnifiedStateModelElements( Rf1HodographQuaternionIndex );
    double Rf2Hodograph = currentUnifiedStateModelElements( Rf2HodographQuaternionIndex );
    double epsilon1Quaternion = currentUnifiedStateModelElements( epsilon1QuaternionIndex );
    double epsilon2Quaternion = currentUnifiedStateModelElements( epsilon2QuaternionIndex );
    double epsilon3Quaternion = currentUnifiedStateModelElements( epsilon3QuaternionIndex );
    double etaQuaternion = currentUnifiedStateModelElements( etaQuaternionIndex );

    // Compute supporting parameters
    double quaterionParameter = std::pow( epsilon3Quaternion, 2 ) + std::pow( etaQuaternion, 2 );
    double sineLambdaParameter = ( 2 * epsilon3Quaternion * etaQuaternion ) / quaterionParameter;
    double cosineLambdaParameter =  ( std::pow( etaQuaternion, 2 ) - std::pow( epsilon3Quaternion, 2 ) ) /
            quaterionParameter;
    double gammaParameter = ( epsilon1Quaternion * epsilon3Quaternion -
                              epsilon2Quaternion * etaQuaternion ) / quaterionParameter;

    double velocityHodographParameter = CHodograph - Rf1Hodograph * sineLambdaParameter +
            Rf2Hodograph * cosineLambdaParameter;
    Eigen::Vector3d rotationalVelocity = Eigen::Vector3d::Zero( );
    rotationalVelocity( 0 ) = accelerationsInRswFrame( 2 ) / velocityHodographParameter;
    rotationalVelocity( 2 ) = std::pow( velocityHodographParameter, 2 ) * CHodograph /
            centralBodyGravitationalParameter;

    Eigen::Vector3d pParameter = Eigen::Vector3d::Zero( );
    pParameter( 0 ) = CHodograph;
    pParameter( 1 ) = Rf2Hodograph;
    pParameter( 2 ) = Rf1Hodograph;
    pParameter = pParameter / velocityHodographParameter;

    // Evaluate USM7 equations
    return computeStateDerivativeForUnifiedStateModelWithQuaternions(
                currentUnifiedStateModelElements, accelerationsInRswFrame, sineLambdaParameter,
                cosineLambdaParameter, gammaParameter, rotationalVelocity, pParameter );
}

//! Function to evaluate the state derivative for the unified state model with quaternions
Eigen::Vector7d computeStateDerivativeForUnifiedStateModelWithQuaternions(
        const Eigen::Vector7d& currentUnifiedStateModelElements,
        const Eigen::Vector6d& currentCartesianState,
        const Eigen::Vector3d& accelerationsInInertialFrame,
        const double centralBodyGravitationalParameter )
{
    return computeStateDerivativeForUnifiedStateModelWithQuaternions(
                currentUnifiedStateModelElements,
                reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrx(
                    currentCartesianState ) * accelerationsInInertialFrame, centralBodyGravitationalParameter );
}


} // namespace propagators

} // namespace tudat

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
        const Eigen::Vector3d& rotationalVelocityVector,
        const Eigen::Vector3d& pParameterVector )
{
    // REMOVE vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv REMOVE
    std::cout << "USM state: " << std::endl << currentUnifiedStateModelElements << std::endl;
    std::cout << "Norm before propagation: " << currentUnifiedStateModelElements.segment( 3, 4 ).norm( ) << std::endl;
    std::cout << "Acceleration: " << std::endl << accelerationsInRswFrame << std::endl;
    std::cout << "Sine lambda: " << sineLambdaParameter << std::endl;
    std::cout << "Cosine lambda: " << cosineLambdaParameter << std::endl;
    std::cout << "Gamma: " << gammaParameter << std::endl;
    std::cout << "Rotational velocity: " << std::endl << rotationalVelocityVector << std::endl;
    std::cout << "P parameter: " << std::endl << pParameterVector << std::endl;
    // REMOVE ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ REMOVE

    // Compute supporting parameters
    Eigen::Matrix3d hodographMatrix = Eigen::Matrix3d::Zero( );
    hodographMatrix( 0, 1 ) = - pParameterVector( 0 );
    hodographMatrix( 1, 0 ) = cosineLambdaParameter;
    hodographMatrix( 1, 1 ) = - ( 1 + pParameterVector( 0 ) ) * sineLambdaParameter;
    hodographMatrix( 1, 2 ) = - gammaParameter * pParameterVector( 1 );
    hodographMatrix( 2, 0 ) = sineLambdaParameter;
    hodographMatrix( 2, 1 ) = ( 1 + pParameterVector( 0 ) ) * cosineLambdaParameter;
    hodographMatrix( 2, 1 ) = gammaParameter * pParameterVector( 2 );

    Eigen::MatrixXd quaternionMatrix = getQuaterionToQuaternionRateMatrix( rotationalVelocityVector );
    // REMOVE vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv REMOVE
    std::cout << "Quaternion augmented matrix: " << std::endl << quaternionMatrix << std::endl;
    // REMOVE ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ REMOVE

    // Evaluate USM7 equations.
    Eigen::Vector7d stateDerivative;
    stateDerivative.segment( 0, 3 ) = hodographMatrix * accelerationsInRswFrame;
    stateDerivative.segment( 3, 4 ) = quaternionMatrix * currentUnifiedStateModelElements.segment( 3, 4 );
    // the 0.5 constant is already accounted for in quaternionMatrix

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
    Eigen::Vector3d rotationalVelocityVector = Eigen::Vector3d::Zero( );
    rotationalVelocityVector( 0 ) = accelerationsInRswFrame( 2 ) / velocityHodographParameter;
    rotationalVelocityVector( 2 ) = std::pow( velocityHodographParameter, 2 ) * CHodograph /
            centralBodyGravitationalParameter;

    Eigen::Vector3d pParameterVector = Eigen::Vector3d::Zero( );
    pParameterVector( 0 ) = CHodograph;
    pParameterVector( 1 ) = Rf2Hodograph;
    pParameterVector( 2 ) = Rf1Hodograph;
    pParameterVector = pParameterVector / velocityHodographParameter;

    // Evaluate USM7 equations
    return computeStateDerivativeForUnifiedStateModelWithQuaternions(
                currentUnifiedStateModelElements, accelerationsInRswFrame, sineLambdaParameter,
                cosineLambdaParameter, gammaParameter, rotationalVelocityVector, pParameterVector );
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

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
#include "Tudat/Astrodynamics/Propagators/nBodyUnifiedStateModelQuaternionsStateDerivative.h"

namespace tudat
{

namespace propagators
{

//! Function to evaluate the state derivative for the unified state model with quaternions
Eigen::Vector7d computeStateDerivativeForUnifiedStateModelQuaternions(
        const Eigen::Vector7d& currentUnifiedStateModelElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double sineLambda,
        const double cosineLambda,
        const double gammaParameter,
        const Eigen::Vector3d& rotationalVelocityVector,
        const Eigen::Vector3d& pAuxiliaryVector )
{
    // Compute supporting parameters
    Eigen::Matrix3d hodographMatrix = Eigen::Matrix3d::Zero( );
    hodographMatrix( 0, 1 ) = - pAuxiliaryVector( 0 );
    hodographMatrix( 1, 0 ) = cosineLambda;
    hodographMatrix( 1, 1 ) = - ( 1.0 + pAuxiliaryVector( 0 ) ) * sineLambda;
    hodographMatrix( 1, 2 ) = - gammaParameter * pAuxiliaryVector( 1 );
    hodographMatrix( 2, 0 ) = sineLambda;
    hodographMatrix( 2, 1 ) = ( 1.0 + pAuxiliaryVector( 0 ) ) * cosineLambda;
    hodographMatrix( 2, 2 ) = gammaParameter * pAuxiliaryVector( 2 );

    Eigen::Matrix4d quaternionMatrix = Eigen::Matrix4d::Zero( );
    // getQuaterionToQuaternionRateMatrix( rotationalVelocityVector.reverse( ) ).reverse( ); // still wrong
    // if the function above is used, remove the 0.5 from line 56 (state derivative)
    quaternionMatrix( 0, 1 ) =   rotationalVelocityVector( 2 );
    quaternionMatrix( 0, 3 ) =   rotationalVelocityVector( 0 );
    quaternionMatrix( 1, 0 ) = - rotationalVelocityVector( 2 );
    quaternionMatrix( 1, 2 ) =   rotationalVelocityVector( 0 );
    quaternionMatrix( 2, 1 ) = - rotationalVelocityVector( 0 );
    quaternionMatrix( 2, 3 ) =   rotationalVelocityVector( 2 );
    quaternionMatrix( 3, 0 ) = - rotationalVelocityVector( 0 );
    quaternionMatrix( 3, 2 ) = - rotationalVelocityVector( 2 );

    // Evaluate USM7 equations.
    Eigen::Vector7d stateDerivative = Eigen::Vector7d::Zero( );
    stateDerivative.segment( 0, 3 ) = hodographMatrix * accelerationsInRswFrame;
    stateDerivative.segment( 3, 4 ) = 0.5 * quaternionMatrix * currentUnifiedStateModelElements.segment( 3, 4 );

    // Give output
    return stateDerivative;
}

//! Function to evaluate the state derivative for the unified state model with quaternions
Eigen::Vector7d computeStateDerivativeForUnifiedStateModelQuaternions(
        const Eigen::Vector7d& currentUnifiedStateModelElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double centralBodyGravitationalParameter )
{
    using namespace orbital_element_conversions;

    // Retrieve USM elements
    double epsilon1Quaternion = currentUnifiedStateModelElements( epsilon1QuaternionIndex );
    double epsilon2Quaternion = currentUnifiedStateModelElements( epsilon2QuaternionIndex );
    double epsilon3Quaternion = currentUnifiedStateModelElements( epsilon3QuaternionIndex );
    double etaQuaternion = currentUnifiedStateModelElements( etaQuaternionIndex );

    // Compute supporting parameters
    double denominator = std::pow( epsilon3Quaternion, 2 ) + std::pow( etaQuaternion, 2 );
    double sineLambda = ( 2 * epsilon3Quaternion * etaQuaternion ) / denominator;
    double cosineLambda =  ( std::pow( etaQuaternion, 2 ) - std::pow( epsilon3Quaternion, 2 ) ) / denominator;
    double gammaParameter = ( epsilon1Quaternion * epsilon3Quaternion -
                              epsilon2Quaternion * etaQuaternion ) / denominator;

    double velocityHodographParameter = currentUnifiedStateModelElements( CHodographQuaternionIndex ) -
            currentUnifiedStateModelElements( Rf1HodographQuaternionIndex ) * sineLambda +
            currentUnifiedStateModelElements( Rf2HodographQuaternionIndex ) * cosineLambda;
    Eigen::Vector3d rotationalVelocityVector = Eigen::Vector3d::Zero( );
    rotationalVelocityVector( 0 ) = accelerationsInRswFrame( 2 ) / velocityHodographParameter;
    rotationalVelocityVector( 2 ) = std::pow( velocityHodographParameter, 2 ) *
            currentUnifiedStateModelElements( CHodographQuaternionIndex ) / centralBodyGravitationalParameter;

    Eigen::Vector3d pAuxiliaryVector = Eigen::Vector3d::Zero( );
    pAuxiliaryVector( 0 ) = currentUnifiedStateModelElements( CHodographQuaternionIndex );
    pAuxiliaryVector( 1 ) = currentUnifiedStateModelElements( Rf2HodographQuaternionIndex );
    pAuxiliaryVector( 2 ) = currentUnifiedStateModelElements( Rf1HodographQuaternionIndex );
    pAuxiliaryVector /= velocityHodographParameter;

    // Evaluate USM7 equations
    return computeStateDerivativeForUnifiedStateModelQuaternions(
                currentUnifiedStateModelElements, accelerationsInRswFrame, sineLambda,
                cosineLambda, gammaParameter, rotationalVelocityVector, pAuxiliaryVector );
}

//! Function to evaluate the state derivative for the unified state model with quaternions
Eigen::Vector7d computeStateDerivativeForUnifiedStateModelQuaternions(
        const Eigen::Vector7d& currentUnifiedStateModelElements,
        const Eigen::Vector6d& currentCartesianState,
        const Eigen::Vector3d& accelerationsInInertialFrame,
        const double centralBodyGravitationalParameter )
{
    return computeStateDerivativeForUnifiedStateModelQuaternions(
                currentUnifiedStateModelElements,
                reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrx(
                    currentCartesianState ) * accelerationsInInertialFrame, centralBodyGravitationalParameter );
}


} // namespace propagators

} // namespace tudat

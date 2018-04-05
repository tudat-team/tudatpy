/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/Propagators/nBodyUnifiedStateModelWithQuaternionsStateDerivative.h"

namespace tudat
{

namespace propagators
{

//! Function to evaluate the state derivative for the unified state model with quaternions
Eigen::VectorXd computeStateDerivativeForUnifiedStateModelWithQuaternions(
        const Eigen::VectorXd& currentUnifiedStateModelElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double sineLambdaParameter,
        const double cosineLambdaParameter,
        const double gammaParameter,
        const Eigen::Vector3d& rotationalVelocity,
        const Eigen::Vector3d& pParameter )
{

    // Compute supporting parameters
    Eigen::Matrix3d hodographMatrix = Eigen::Matrix3d::Zero( );
    hodographMatrix( 1, 2 ) = - pParameter( 0 );
    hodographMatrix( 2, 1 ) = cosineLambdaParameter;
    hodographMatrix( 2, 2 ) = - ( 1 + pParameter( 0 ) ) * sineLambdaParameter;
    hodographMatrix( 2, 3 ) = - gammaParameter * pParameter( 1 );
    hodographMatrix( 3, 1 ) = sineLambdaParameter;
    hodographMatrix( 3, 2 ) = ( 1 + pParameter( 0 ) ) * cosineLambdaParameter;
    hodographMatrix( 3, 3 ) = gammaParameter * pParameter( 2 );

    Eigen::MatrixXd quaternionMatrix = Eigen::MatrixXd::Zero( 4, 4 );
    quaternionMatrix( 1, 2 ) = rotationalVelocity( 2 );
    quaternionMatrix( 1, 3 ) = - rotationalVelocity( 1 );
    quaternionMatrix( 1, 4 ) = rotationalVelocity( 0 );
    quaternionMatrix( 2, 1 ) = - rotationalVelocity( 2 );
    quaternionMatrix( 2, 3 ) = rotationalVelocity( 0 );
    quaternionMatrix( 2, 4 ) = rotationalVelocity( 1 );
    quaternionMatrix( 3, 1 ) = rotationalVelocity( 1 );
    quaternionMatrix( 3, 2 ) = - rotationalVelocity( 0 );
    quaternionMatrix( 3, 4 ) = rotationalVelocity( 2 );
    quaternionMatrix( 4, 1 ) = - rotationalVelocity( 0 );
    quaternionMatrix( 4, 2 ) = - rotationalVelocity( 1 );
    quaternionMatrix( 4, 3 ) = - rotationalVelocity( 2 );

    Eigen::VectorXd stateDerivative;

    // Evaluate USM7 equations.
    stateDerivative.segment( 0, 3 ) = hodographMatrix * accelerationsInRswFrame;
    stateDerivative.segment( 3, 3 ) = quaternionMatrix * currentUnifiedStateModelElements.segment( 3, 3 );

    return stateDerivative;

}

//! Function to evaluate the state derivative for the unified state model with quaternions
Eigen::VectorXd computeStateDerivativeForUnifiedStateModelWithQuaternions(
        const Eigen::VectorXd& currentUnifiedStateModelElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double centralBodyGravitationalParameter )
{

    // Retrieve USM elements
    double CHodograph = currentUnifiedStateModelElements( 0 );
    double Rf1Hodograph = currentUnifiedStateModelElements( 1 );
    double Rf2Hodograph = currentUnifiedStateModelElements( 2 );
    double epsilon1Quaternion = currentUnifiedStateModelElements( 3 );
    double epsilon2Quaternion = currentUnifiedStateModelElements( 4 );
    double epsilon3Quaternion = currentUnifiedStateModelElements( 5 );
    double etaQuaternion = currentUnifiedStateModelElements( 6 );

    // Compute supporting parameters
    double quaterionParameter = std::pow( epsilon3Quaternion, 2) + std::pow( etaQuaternion, 2 );
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

    Eigen::Vector3d pParameter;
    pParameter( 0 ) = CHodograph;
    pParameter( 1 ) = Rf1Hodograph;
    pParameter( 2 ) = Rf2Hodograph;
    pParameter = pParameter / velocityHodographParameter;

    // Evaluate USM7 equations
    return computeStateDerivativeForUnifiedStateModelWithQuaternions(
                currentUnifiedStateModelElements, accelerationsInRswFrame, sineLambdaParameter,
                cosineLambdaParameter, gammaParameter, rotationalVelocity, pParameter );
}

//! Function to evaluate the state derivative for the unified state model with quaternions
Eigen::VectorXd computeStateDerivativeForUnifiedStateModelWithQuaternions(
        const Eigen::VectorXd& currentUnifiedStateModelElements,
        const Eigen::Vector6d& currentCartesianState,
        const Eigen::Vector3d& accelerationsInInertialFrame,
        const double centralBodyGravitationalParameter )
{

    return  computeStateDerivativeForUnifiedStateModelWithQuaternions(
                currentUnifiedStateModelElements,
                reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrx(
                    currentCartesianState ) * accelerationsInInertialFrame, centralBodyGravitationalParameter );
}


} // namespace propagators

} // namespace tudat

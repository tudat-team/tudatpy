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
#include "Tudat/Astrodynamics/Propagators/nBodyUnifiedStateModelWithModifiedRodriguesParametersStateDerivative.h"
//#include "Tudat/Astrodynamics/Propagators/rotationalMotionStateDerivative.h"

namespace tudat
{

namespace propagators
{

//! Function to evaluate the state derivative for the unified state model with modified rodrigues parameters
Eigen::Vector7d computeStateDerivativeForUnifiedStateModelWithModifiedRodriguesParameters(
        const Eigen::Vector7d& currentUnifiedStateModelElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double sineLambdaParameter,
        const double cosineLambdaParameter,
        const double gammaParameter,
        const Eigen::Vector3d& rotationalVelocityVector,
        const Eigen::Vector3d& pParameterVector )
{
    // Compute supporting parameters
    Eigen::Matrix3d hodographMatrix = Eigen::Matrix3d::Zero( );
    hodographMatrix( 0, 1 ) = - pParameterVector( 0 );
    hodographMatrix( 1, 0 ) = cosineLambdaParameter;
    hodographMatrix( 1, 1 ) = - ( 1.0 + pParameterVector( 0 ) ) * sineLambdaParameter;
    hodographMatrix( 1, 2 ) = - gammaParameter * pParameterVector( 1 );
    hodographMatrix( 2, 0 ) = sineLambdaParameter;
    hodographMatrix( 2, 1 ) = ( 1.0 + pParameterVector( 0 ) ) * cosineLambdaParameter;
    hodographMatrix( 2, 2 ) = gammaParameter * pParameterVector( 2 );

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

    // Evaluate USM6 equations.
    Eigen::Vector7d stateDerivative;
    stateDerivative.segment( 0, 3 ) = hodographMatrix * accelerationsInRswFrame;
    stateDerivative.segment( 3, 4 ) = 0.5 * quaternionMatrix * currentUnifiedStateModelElements.segment( 3, 4 );

    // Give output
    return stateDerivative;
}

//! Function to evaluate the state derivative for the unified state model with modified rodrigues parameters
Eigen::Vector7d computeStateDerivativeForUnifiedStateModelWithModifiedRodriguesParameters(
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

    // Evaluate USM6 equations
    return computeStateDerivativeForUnifiedStateModelWithModifiedRodriguesParameters(
                currentUnifiedStateModelElements, accelerationsInRswFrame, sineLambdaParameter,
                cosineLambdaParameter, gammaParameter, rotationalVelocityVector, pParameterVector );
}

//! Function to evaluate the state derivative for the unified state model with modified rodrigues parameters
Eigen::Vector7d computeStateDerivativeForUnifiedStateModelWithModifiedRodriguesParameters(
        const Eigen::Vector7d& currentUnifiedStateModelElements,
        const Eigen::Vector6d& currentCartesianState,
        const Eigen::Vector3d& accelerationsInInertialFrame,
        const double centralBodyGravitationalParameter )
{
    return computeStateDerivativeForUnifiedStateModelWithModifiedRodriguesParameters(
                currentUnifiedStateModelElements,
                reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrx(
                    currentCartesianState ) * accelerationsInInertialFrame, centralBodyGravitationalParameter );
}


} // namespace propagators

} // namespace tudat

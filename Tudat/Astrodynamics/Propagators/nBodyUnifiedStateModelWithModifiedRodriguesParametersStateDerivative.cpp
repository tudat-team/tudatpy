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

namespace tudat
{

namespace propagators
{

//! Function to evaluate the state derivative for the unified state model with modified rodrigues parameters
Eigen::Vector7d computeStateDerivativeForUnifiedStateModelWithModifiedRodriguesParameters(
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

    // Compute kinematic equation, i.e., derivative of modified Rodrigues parameters (also valid for SMRP)
    Eigen::Vector3d modifiedRodriguesParametersVector = currentUnifiedStateModelElements.segment( 3, 3 );
    Eigen::Matrix3d skewModifiedRodriguesParametersVector = linear_algebra::getCrossProductMatrix( modifiedRodriguesParametersVector );
    Eigen::Vector3d modifiedRodriguesParametersDerivative = 0.5 *
            ( 0.5 * ( 1.0 - std::pow( modifiedRodriguesParametersVector.norm( ), 2 ) ) * Eigen::Matrix3d::Identity( ) +
              skewModifiedRodriguesParametersVector + modifiedRodriguesParametersVector * modifiedRodriguesParametersVector.transpose( ) ) *
            rotationalVelocityVector;

    // Evaluate USM6 equations.
    Eigen::Vector7d stateDerivative = Eigen::Vector7d::Zero( );
    stateDerivative.segment( 0, 3 ) = hodographMatrix * accelerationsInRswFrame;
    stateDerivative.segment( 3, 3 ) = modifiedRodriguesParametersDerivative;

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

    // Compute auxiliary parameters
    Eigen::Vector3d modifiedRodriguesParametersVector = currentUnifiedStateModelElements.segment( sigma1ModifiedRodriguesParameterIndex, 3 );
    double modifiedRodriguesParametersMagnitude = modifiedRodriguesParametersVector.norm( );
    // magnitude of modified rodrigues parameters, also called sigma

    // Precompute often used variables
    // Note the for SMRP some variable names do not match their definitions
    double modifiedRodriguesParametersMagnitudeSquared =
            modifiedRodriguesParametersMagnitude * modifiedRodriguesParametersMagnitude;
    double oneMinusModifiedRodriguesParametersMagnitudeSquared =
            currentUnifiedStateModelElements( shadowModifiedRodriguesParameterFlagIndex ) ?
                ( modifiedRodriguesParametersMagnitudeSquared - 1.0 ) : // inverse definition for SMRP
                ( 1.0 - modifiedRodriguesParametersMagnitudeSquared );
    double oneMinusModifiedRodriguesParametersMagnitudeSquaredSquared =
            std::pow( oneMinusModifiedRodriguesParametersMagnitudeSquared, 2 );
    double sigma1ModifiedRodriguesParametersSquared = std::pow( modifiedRodriguesParametersVector( 0 ), 2 );
    double sigma2ModifiedRodriguesParametersSquared = std::pow( modifiedRodriguesParametersVector( 1 ), 2 );
    double sigma3ModifiedRodriguesParametersSquared = std::pow( modifiedRodriguesParametersVector( 2 ), 2 );

    // Compute supporting parameters
    double denominator = 4.0 * sigma3ModifiedRodriguesParametersSquared +
            oneMinusModifiedRodriguesParametersMagnitudeSquaredSquared; // denominator is never null
    double sineLambda = 4.0 * modifiedRodriguesParametersVector( 2 ) *
            ( 1.0 - modifiedRodriguesParametersMagnitudeSquared ) / denominator;
    double cosineLambda = ( oneMinusModifiedRodriguesParametersMagnitudeSquaredSquared -
                            4.0 * sigma3ModifiedRodriguesParametersSquared ) / denominator;
    double gammaParameter = 2.0 * ( modifiedRodriguesParametersVector( 1 ) * ( modifiedRodriguesParametersMagnitudeSquared - 1.0 ) +
                                    2.0 * modifiedRodriguesParametersVector( 0 ) * modifiedRodriguesParametersVector( 2 ) ) /
            ( std::pow( sigma1ModifiedRodriguesParametersSquared, 2 ) + std::pow( sigma2ModifiedRodriguesParametersSquared, 2 ) +
              std::pow( sigma3ModifiedRodriguesParametersSquared, 2 ) + 2.0 *
              ( sigma1ModifiedRodriguesParametersSquared * sigma2ModifiedRodriguesParametersSquared +
                sigma1ModifiedRodriguesParametersSquared * sigma3ModifiedRodriguesParametersSquared +
                sigma2ModifiedRodriguesParametersSquared * sigma3ModifiedRodriguesParametersSquared -
                modifiedRodriguesParametersMagnitudeSquared + 2.0 * sigma3ModifiedRodriguesParametersSquared ) + 1.0 );

    double velocityHodographParameter = currentUnifiedStateModelElements( CHodographModifiedRodriguesParameterIndex ) -
            currentUnifiedStateModelElements( Rf1HodographModifiedRodriguesParameterIndex ) * sineLambda +
            currentUnifiedStateModelElements( Rf2HodographModifiedRodriguesParameterIndex ) * cosineLambda;
    Eigen::Vector3d rotationalVelocityVector = Eigen::Vector3d::Zero( );
    rotationalVelocityVector( 0 ) = accelerationsInRswFrame( 2 ) / velocityHodographParameter;
    rotationalVelocityVector( 2 ) = std::pow( velocityHodographParameter, 2 ) *
            currentUnifiedStateModelElements( CHodographModifiedRodriguesParameterIndex ) / centralBodyGravitationalParameter;

    Eigen::Vector3d pAuxiliaryVector = Eigen::Vector3d::Zero( );
    pAuxiliaryVector( 0 ) = currentUnifiedStateModelElements( CHodographModifiedRodriguesParameterIndex );
    pAuxiliaryVector( 1 ) = currentUnifiedStateModelElements( Rf2HodographModifiedRodriguesParameterIndex );
    pAuxiliaryVector( 2 ) = currentUnifiedStateModelElements( Rf1HodographModifiedRodriguesParameterIndex );
    pAuxiliaryVector /= velocityHodographParameter;

    // Evaluate USM6 equations
    return computeStateDerivativeForUnifiedStateModelWithModifiedRodriguesParameters(
                currentUnifiedStateModelElements, accelerationsInRswFrame, sineLambda,
                cosineLambda, gammaParameter, rotationalVelocityVector, pAuxiliaryVector );
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

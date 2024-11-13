/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Vittaldev, V. (2010). The unified state model: Derivation and application in astro
 *          and navigation. Master thesis, Delft University of Technology.
 */

#include "tudat/astro/basic_astro/attitudeElementConversions.h"

#include "tudat/astro/propagators/nBodyUnifiedStateModelExponentialMapStateDerivative.h"
#include "tudat/astro/propagators/rotationalMotionExponentialMapStateDerivative.h"

namespace tudat
{

namespace propagators
{

//! Function to evaluate the state derivative for the unified state model with exponential map
Eigen::Vector7d computeStateDerivativeForUnifiedStateModelExponentialMap(
        const Eigen::Vector7d& currentUnifiedStateModelElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double sineLambda,
        const double cosineLambda,
        const double gammaParameter,
        const Eigen::Vector3d& rotationalVelocityVector,
        const Eigen::Vector3d& pAuxiliaryVector )
{
    // Compute matrix for dynamic equation
    Eigen::Matrix3d hodographMatrix = Eigen::Matrix3d::Zero( );
    hodographMatrix( 0, 1 ) = - pAuxiliaryVector( 0 );
    hodographMatrix( 1, 0 ) = cosineLambda;
    hodographMatrix( 1, 1 ) = - ( 1.0 + pAuxiliaryVector( 0 ) ) * sineLambda;
    hodographMatrix( 1, 2 ) = - gammaParameter * pAuxiliaryVector( 1 );
    hodographMatrix( 2, 0 ) = sineLambda;
    hodographMatrix( 2, 1 ) = ( 1.0 + pAuxiliaryVector( 0 ) ) * cosineLambda;
    hodographMatrix( 2, 2 ) = gammaParameter * pAuxiliaryVector( 2 );

    // Evaluate USMEM equations.
    Eigen::Vector7d stateDerivative;
    stateDerivative.segment( 0, 3 ) = hodographMatrix * accelerationsInRswFrame;
    stateDerivative.segment( 3, 4 ) = calculateExponentialMapDerivative( currentUnifiedStateModelElements.segment( 3, 4 ),
                                                                         rotationalVelocityVector );

    // Give output
    return stateDerivative;
}

//! Function to evaluate the state derivative for the unified state model with exponential map
Eigen::Vector7d computeStateDerivativeForUnifiedStateModelExponentialMap(
        const Eigen::Vector7d& currentUnifiedStateModelElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double centralBodyGravitationalParameter )
{
    using namespace orbital_element_conversions;

    // Convert exponential map to quaternions
    Eigen::Vector4d quaternionElements = convertExponentialMapToQuaternionElements(
                currentUnifiedStateModelElements.segment( e1USMEMIndex, 4 ) );
    double etaQuaternionParameter = quaternionElements( etaQuaternionIndex );
    Eigen::Vector3d epsilonQuaternionVector = quaternionElements.segment( epsilon1QuaternionIndex, 3 );

    // Compute supporting parameters
    double denominator = std::pow( epsilonQuaternionVector( 2 ), 2 ) + std::pow( etaQuaternionParameter, 2 );
    double sineLambda = ( 2 * epsilonQuaternionVector( 2 ) * etaQuaternionParameter ) / denominator;
    double cosineLambda =  ( std::pow( etaQuaternionParameter, 2 ) - std::pow( epsilonQuaternionVector( 2 ), 2 ) ) / denominator;
    double gammaParameter = ( epsilonQuaternionVector( 0 ) * epsilonQuaternionVector( 2 ) -
                              epsilonQuaternionVector( 1 ) * etaQuaternionParameter ) / denominator;

    double velocityHodographParameter = currentUnifiedStateModelElements( CHodographUSMEMIndex ) -
            currentUnifiedStateModelElements( Rf1HodographUSMEMIndex ) * sineLambda +
            currentUnifiedStateModelElements( Rf2HodographUSMEMIndex ) * cosineLambda;
    Eigen::Vector3d rotationalVelocityVector = Eigen::Vector3d::Zero( );
    rotationalVelocityVector( 0 ) = accelerationsInRswFrame( 2 ) / velocityHodographParameter;
    rotationalVelocityVector( 2 ) = std::pow( velocityHodographParameter, 2 ) *
            currentUnifiedStateModelElements( CHodographUSMEMIndex ) / centralBodyGravitationalParameter;

    Eigen::Vector3d pAuxiliaryVector = Eigen::Vector3d::Zero( );
    pAuxiliaryVector( 0 ) = currentUnifiedStateModelElements( CHodographUSMEMIndex );
    pAuxiliaryVector( 1 ) = currentUnifiedStateModelElements( Rf2HodographUSMEMIndex );
    pAuxiliaryVector( 2 ) = currentUnifiedStateModelElements( Rf1HodographUSMEMIndex );
    pAuxiliaryVector /= velocityHodographParameter;

    // Evaluate USMEM equations
    return computeStateDerivativeForUnifiedStateModelExponentialMap(
                currentUnifiedStateModelElements, accelerationsInRswFrame, sineLambda,
                cosineLambda, gammaParameter, rotationalVelocityVector, pAuxiliaryVector );
}

//! Function to evaluate the state derivative for the unified state model with exponential map
Eigen::Vector7d computeStateDerivativeForUnifiedStateModelExponentialMap(
        const Eigen::Vector7d& currentUnifiedStateModelElements,
        const Eigen::Vector6d& currentCartesianState,
        const Eigen::Vector3d& accelerationsInInertialFrame,
        const double centralBodyGravitationalParameter )
{
    return computeStateDerivativeForUnifiedStateModelExponentialMap(
                currentUnifiedStateModelElements,
                reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrix(
                    currentCartesianState ) * accelerationsInInertialFrame, centralBodyGravitationalParameter );
}

//template class NBodyUnifiedStateModelExponentialMapStateDerivative< double, double >;


} // namespace propagators

} // namespace tudat

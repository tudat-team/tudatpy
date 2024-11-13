/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/basic_astro/stateVectorIndices.h"
#include "tudat/astro/propagators/nBodyUnifiedStateModelQuaternionsStateDerivative.h"
#include "tudat/astro/propagators/rotationalMotionQuaternionsStateDerivative.h"

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

    // Evaluate USM7 equations.
    Eigen::Vector7d stateDerivative;
    stateDerivative.segment( 0, 3 ) = hodographMatrix * accelerationsInRswFrame;
    stateDerivative.segment( 3, 4 ) = calculateQuaternionDerivative( currentUnifiedStateModelElements.segment( 3, 4 ),
                                                                      rotationalVelocityVector );

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
    double etaQuaternionParameter = currentUnifiedStateModelElements( etaUSM7Index );
    double epsilon1QuaternionParameter = currentUnifiedStateModelElements( epsilon1USM7Index );
    double epsilon2QuaternionParameter = currentUnifiedStateModelElements( epsilon2USM7Index );
    double epsilon3QuaternionParameter = currentUnifiedStateModelElements( epsilon3USM7Index );

    // Compute supporting parameters
    double denominator = std::pow( epsilon3QuaternionParameter, 2 ) + std::pow( etaQuaternionParameter, 2 );
    double sineLambda = ( 2 * epsilon3QuaternionParameter * etaQuaternionParameter ) / denominator;
    double cosineLambda =  ( std::pow( etaQuaternionParameter, 2 ) - std::pow( epsilon3QuaternionParameter, 2 ) ) / denominator;
    double gammaParameter = ( epsilon1QuaternionParameter * epsilon3QuaternionParameter -
                              epsilon2QuaternionParameter * etaQuaternionParameter ) / denominator;

    double velocityHodographParameter = currentUnifiedStateModelElements( CHodographUSM7Index ) -
            currentUnifiedStateModelElements( Rf1HodographUSM7Index ) * sineLambda +
            currentUnifiedStateModelElements( Rf2HodographUSM7Index ) * cosineLambda;
    Eigen::Vector3d rotationalVelocityVector = Eigen::Vector3d::Zero( );
    rotationalVelocityVector( 0 ) = accelerationsInRswFrame( 2 ) / velocityHodographParameter;
    rotationalVelocityVector( 2 ) = std::pow( velocityHodographParameter, 2 ) *
            currentUnifiedStateModelElements( CHodographUSM7Index ) / centralBodyGravitationalParameter;

    Eigen::Vector3d pAuxiliaryVector = Eigen::Vector3d::Zero( );
    pAuxiliaryVector( 0 ) = currentUnifiedStateModelElements( CHodographUSM7Index );
    pAuxiliaryVector( 1 ) = currentUnifiedStateModelElements( Rf2HodographUSM7Index );
    pAuxiliaryVector( 2 ) = currentUnifiedStateModelElements( Rf1HodographUSM7Index );
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
                reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrix(
                    currentCartesianState ) * accelerationsInInertialFrame, centralBodyGravitationalParameter );
}

//template class NBodyUnifiedStateModelQuaternionsStateDerivative< double, double >;


} // namespace propagators

} // namespace tudat

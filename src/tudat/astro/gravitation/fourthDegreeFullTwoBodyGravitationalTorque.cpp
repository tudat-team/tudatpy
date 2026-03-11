/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/gravitation/fourthDegreeFullTwoBodyGravitationalTorque.h"

#include <tuple>
#include <vector>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/gravitation/fullTwoBodySphericalHarmonicAcceleration.h"
#include "tudat/astro/gravitation/fullTwoBodySphericalHarmonicTorque.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"

namespace tudat
{

namespace gravitation
{

namespace
{

std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >
getFullDegreeTwoCouplingCombinations( )
{
    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse;
    for( unsigned int m = 0; m <= 2; m++ )
    {
        coefficientCombinationsToUse.push_back( std::make_tuple( 2, m, 0, 0 ) );
    }
    for( unsigned int m = 0; m <= 2; m++ )
    {
        for( unsigned int k = 0; k <= 2; k++ )
        {
            coefficientCombinationsToUse.push_back( std::make_tuple( 2, m, 2, k ) );
        }
    }
    return coefficientCombinationsToUse;
}

}  // namespace

Eigen::Vector3d calculateFourthDegreeFullTwoBodyGravitationalTorque(
        const Eigen::Vector3d& relativePositionOfBodyExertingTorqueInBodyFixedFrame,
        const double massOfBodyExertingTorque,
        const Eigen::Matrix3d& inertiaTensorOfBodyUndergoingTorque,
        const Eigen::Matrix3d& inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque )
{
    const double relativeDistance = relativePositionOfBodyExertingTorqueInBodyFixedFrame.norm( );
    if( relativeDistance <= 0.0 )
    {
        throw std::runtime_error(
                "Error when computing fourth-degree full two-body gravitational torque: relative distance is zero." );
    }

    const double gravitationalParameterOfBodyExertingTorque =
            physical_constants::GRAVITATIONAL_CONSTANT * massOfBodyExertingTorque;
    const double referenceRadius = 1.0;

    const std::tuple< Eigen::MatrixXd, Eigen::MatrixXd, double > body1DegreeTwoCoefficients =
            getDegreeTwoSphericalHarmonicCoefficients(
                    inertiaTensorOfBodyUndergoingTorque,
                    gravitationalParameterOfBodyExertingTorque,
                    referenceRadius,
                    2,
                    true );
    const std::tuple< Eigen::MatrixXd, Eigen::MatrixXd, double > body2DegreeTwoCoefficients =
            getDegreeTwoSphericalHarmonicCoefficients(
                    inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque,
                    gravitationalParameterOfBodyExertingTorque,
                    referenceRadius,
                    2,
                    true );

    const Eigen::MatrixXd cosineCoefficientsOfBody1 = std::get< 0 >( body1DegreeTwoCoefficients );
    const Eigen::MatrixXd sineCoefficientsOfBody1 = std::get< 1 >( body1DegreeTwoCoefficients );
    const Eigen::MatrixXd cosineCoefficientsOfBody2 = std::get< 0 >( body2DegreeTwoCoefficients );
    const Eigen::MatrixXd sineCoefficientsOfBody2 = std::get< 1 >( body2DegreeTwoCoefficients );

    const Eigen::Vector3d positionOfBodyUndergoingTorque = Eigen::Vector3d::Zero( );
    const Eigen::Vector3d positionOfBodyExertingTorque = relativePositionOfBodyExertingTorqueInBodyFixedFrame;
    const Eigen::Quaterniond identityRotation = Eigen::Quaterniond::Identity( );

    // Evaluate Eq. (11) through the equivalent degree-2 mutual interaction terms.
    const std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > degreeTwoMutualAcceleration =
            std::make_shared< FullTwoBodySphericalHarmonicAcceleration >(
                    [ = ]( ) { return positionOfBodyUndergoingTorque; },
                    [ = ]( ) { return positionOfBodyExertingTorque; },
                    [ = ]( ) { return gravitationalParameterOfBodyExertingTorque; },
                    referenceRadius,
                    referenceRadius,
                    [ = ]( ) { return cosineCoefficientsOfBody1; },
                    [ = ]( ) { return sineCoefficientsOfBody1; },
                    [ = ]( ) { return cosineCoefficientsOfBody2; },
                    [ = ]( ) { return sineCoefficientsOfBody2; },
                    getFullDegreeTwoCouplingCombinations( ),
                    [ = ]( ) { return identityRotation; },
                    [ = ]( ) { return identityRotation; },
                    false,
                    true );

    FullTwoBodySphericalHarmonicTorque degreeTwoMutualTorque( degreeTwoMutualAcceleration, true );
    degreeTwoMutualTorque.updateMembers( 0.0 );

    return massOfBodyExertingTorque * degreeTwoMutualTorque.getTorque( );
}

FourthDegreeFullTwoBodyGravitationalTorqueModel::FourthDegreeFullTwoBodyGravitationalTorqueModel(
        const std::function< Eigen::Vector3d( ) > positionOfBodyUndergoingTorqueFunction,
        const std::function< Eigen::Vector3d( ) > positionOfBodyExertingTorqueFunction,
        const std::function< double( ) > massOfBodyExertingTorqueFunction,
        const std::function< Eigen::Matrix3d( ) > inertiaTensorOfBodyUndergoingTorqueFunction,
        const std::function< Eigen::Matrix3d( ) > inertiaTensorOfBodyExertingTorqueFunction,
        const std::function< Eigen::Quaterniond( ) > rotationToBodyFixedFrameOfBodyUndergoingTorqueFunction,
        const std::function< Eigen::Quaterniond( ) > rotationToBodyFixedFrameOfBodyExertingTorqueFunction ):
    positionOfBodyUndergoingTorqueFunction_( positionOfBodyUndergoingTorqueFunction ),
    positionOfBodyExertingTorqueFunction_( positionOfBodyExertingTorqueFunction ),
    massOfBodyExertingTorqueFunction_( massOfBodyExertingTorqueFunction ),
    inertiaTensorOfBodyUndergoingTorqueFunction_( inertiaTensorOfBodyUndergoingTorqueFunction ),
    inertiaTensorOfBodyExertingTorqueFunction_( inertiaTensorOfBodyExertingTorqueFunction ),
    rotationToBodyFixedFrameOfBodyUndergoingTorqueFunction_( rotationToBodyFixedFrameOfBodyUndergoingTorqueFunction ),
    rotationToBodyFixedFrameOfBodyExertingTorqueFunction_( rotationToBodyFixedFrameOfBodyExertingTorqueFunction )
{
}

void FourthDegreeFullTwoBodyGravitationalTorqueModel::updateMembers( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        currentRotationToBodyFixedFrameOfBodyUndergoingTorque_ = rotationToBodyFixedFrameOfBodyUndergoingTorqueFunction_( );
        currentRotationToBodyFixedFrameOfBodyExertingTorque_ = rotationToBodyFixedFrameOfBodyExertingTorqueFunction_( );
        currentRelativePositionInInertialFrame_ =
                positionOfBodyExertingTorqueFunction_( ) - positionOfBodyUndergoingTorqueFunction_( );
        currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_ =
                currentRotationToBodyFixedFrameOfBodyUndergoingTorque_ * currentRelativePositionInInertialFrame_;

        currentMassOfBodyExertingTorque_ = massOfBodyExertingTorqueFunction_( );
        currentInertiaTensorOfBodyUndergoingTorque_ = inertiaTensorOfBodyUndergoingTorqueFunction_( );
        currentInertiaTensorOfBodyExertingTorque_ = inertiaTensorOfBodyExertingTorqueFunction_( );

        const Eigen::Matrix3d rotationFromBody2ToBody1 =
                currentRotationToBodyFixedFrameOfBodyUndergoingTorque_.toRotationMatrix( ) *
                currentRotationToBodyFixedFrameOfBodyExertingTorque_.toRotationMatrix( ).transpose( );
        currentInertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque_ =
                rotationFromBody2ToBody1 * currentInertiaTensorOfBodyExertingTorque_ * rotationFromBody2ToBody1.transpose( );

        currentTorque_ = calculateFourthDegreeFullTwoBodyGravitationalTorque(
                currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_,
                currentMassOfBodyExertingTorque_,
                currentInertiaTensorOfBodyUndergoingTorque_,
                currentInertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque_ );

        currentTime_ = currentTime;
    }
}

}  // namespace gravitation

}  // namespace tudat

/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/RotationalDynamicsPartials/sphericalHarmonicGravitationalTorquePartial.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"


namespace tudat
{

namespace acceleration_partials
{

//! Function to compute partial derivative of spherical harmonic torque w.r.t. quaternion elements
Eigen::Matrix< double, 3, 4 > getPartialDerivativeOfSphericalHarmonicGravitationalTorqueWrtQuaternion(
        const Eigen::Matrix3d& bodyFixedRelativePositionCrossProductMatrix,
        const Eigen::Matrix3d& bodyFixedPotentialGradientPositionPartial,
        const Eigen::Matrix3d& bodyFixedPotentialGradientCrossProductMatrix,
        const Eigen::Vector3d& inertialRelativePosition,
        const std::vector< Eigen::Matrix3d > derivativeOfRotationMatrixWrtQuaternions )
{
    Eigen::Matrix< double, 3, 4 > partialDerivative = Eigen::Matrix< double, 3, 4 >::Zero( );
    for( unsigned int i = 0; i < derivativeOfRotationMatrixWrtQuaternions.size( ); i++ )
    {
        partialDerivative.block( 0, i, 3, 1 ) =
                derivativeOfRotationMatrixWrtQuaternions.at( i ).transpose( ) * inertialRelativePosition;
    }
    return ( bodyFixedRelativePositionCrossProductMatrix * bodyFixedPotentialGradientPositionPartial -
             bodyFixedPotentialGradientCrossProductMatrix ) * partialDerivative;
}

//! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
SphericalHarmonicGravitationalTorquePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionPair;
    partialFunctionPair = std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );

    if( !estimatable_parameters::isParameterRotationMatrixProperty( parameter->getParameterName( ).first ) )
    {
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > accelerationPartialFunction =
                accelerationPartial_->getParameterPartialFunction( parameter );
        if( accelerationPartialFunction.second > 0 )
        {
            partialFunctionPair = std::make_pair(
                        std::bind( &SphericalHarmonicGravitationalTorquePartial::getParameterPartialFromAccelerationPartialFunction,
                                     this, std::placeholders::_1, accelerationPartialFunction ), accelerationPartialFunction.second );
        }
    }

    return partialFunctionPair;
}

//! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > SphericalHarmonicGravitationalTorquePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionPair;
    partialFunctionPair = std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );

    if( !estimatable_parameters::isParameterRotationMatrixProperty( parameter->getParameterName( ).first ) )
    {
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > accelerationPartialFunction =
                accelerationPartial_->getParameterPartialFunction( parameter );
        if( accelerationPartialFunction.second > 0 )
        {
            partialFunctionPair = std::make_pair(
                        std::bind( &SphericalHarmonicGravitationalTorquePartial::getParameterPartialFromAccelerationPartialFunction,
                                     this, std::placeholders::_1, accelerationPartialFunction ), accelerationPartialFunction.second );
        }
    }

    return partialFunctionPair;
}

//! Function for calculating the partial of the torque w.r.t. a non-rotational integrated state
void SphericalHarmonicGravitationalTorquePartial::wrtNonRotationalStateOfAdditionalBody(
        Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const std::pair< std::string, std::string >& stateReferencePoint,
        const propagators::IntegratedStateType integratedStateType )
{
    if( ( stateReferencePoint.first == bodyExertingTorque_ ||
          stateReferencePoint.first == bodyUndergoingTorque_ ) && integratedStateType == propagators::translational_state )
    {
        partialMatrix.block( 0, 0, 3, 3 ) +=
                ( ( stateReferencePoint.first == bodyExertingTorque_ ) ? 1.0 : -1.0 ) *
                currentBodyFixedRelativePositionCrossProductMatrix_ * currentRotationToBodyFixedFrame_ *
                                accelerationPartial_->getCurrentPartialWrtPosition( ) -
                                currentBodyFixedPotentialGradientCrossProductMatrix_ * currentRotationToBodyFixedFrame_;
    }
}

//! Update partial model to current time
void SphericalHarmonicGravitationalTorquePartial::update( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        torqueModel_->updateMembers( currentTime );
        accelerationPartial_->update( currentTime );

        currentRotationToBodyFixedFrame_ =
                torqueModel_->getSphericalHarmonicAcceleration( )->getCurrentRotationToIntegrationFrameMatrix( ).transpose( );

        currentBodyFixedRelativePosition_ = torqueModel_->getSphericalHarmonicAcceleration( )->getCurrentRelativePosition( );
        currentBodyFixedRelativePositionCrossProductMatrix_ = linear_algebra::getCrossProductMatrix(
                    currentBodyFixedRelativePosition_ );
        currentBodyFixedPotentialGradient_ = torqueModel_->getSphericalHarmonicAcceleration( )->getAccelerationInBodyFixedFrame( );
        currentBodyFixedPotentialGradientCrossProductMatrix_ = linear_algebra::getCrossProductMatrix(
                    currentBodyFixedPotentialGradient_ );

        currentQuaternionVector_ = linear_algebra::convertQuaternionToVectorFormat(
                    Eigen::Quaterniond( currentRotationToBodyFixedFrame_.transpose( ) ) );
        linear_algebra::computePartialDerivativeOfRotationMatrixWrtQuaternion(
                    currentQuaternionVector_,  currentRotationMatrixDerivativesWrtQuaternion_ );

        currentPartialDerivativeWrtQuaternion_ =
                getPartialDerivativeOfSphericalHarmonicGravitationalTorqueWrtQuaternion(
                    currentBodyFixedRelativePositionCrossProductMatrix_,
                    accelerationPartial_->getCurrentBodyFixedPartialWrtPosition( ),
                    currentBodyFixedPotentialGradientCrossProductMatrix_,
                    torqueModel_->getSphericalHarmonicAcceleration( )->getCurrentInertialRelativePosition( ),
                    currentRotationMatrixDerivativesWrtQuaternion_ );

        currentParameterPartialPreMultiplier_ = currentBodyFixedRelativePositionCrossProductMatrix_ *
                currentRotationToBodyFixedFrame_;
    }
}

void SphericalHarmonicGravitationalTorquePartial::getParameterPartialFromAccelerationPartialFunction(
        Eigen::MatrixXd& partialMatrix,
        const std::pair< std::function< void( Eigen::MatrixXd& ) >, int >& accelerationPartialFunction )
{
    Eigen::MatrixXd accelerationPartialsMatrix = Eigen::MatrixXd( 3, accelerationPartialFunction.second );
    accelerationPartialFunction.first( accelerationPartialsMatrix );

    partialMatrix += currentParameterPartialPreMultiplier_ * accelerationPartialsMatrix;
}

}

}

/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/rotational_dynamics_partials/fullTwoBodySphericalHarmonicGravitationalTorquePartial.h"

#include "tudat/astro/orbit_determination/acceleration_partials/numericalAccelerationPartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"

namespace tudat
{

namespace acceleration_partials
{

FullTwoBodySphericalHarmonicGravitationalTorquePartial::FullTwoBodySphericalHarmonicGravitationalTorquePartial(
        const std::shared_ptr< gravitation::FullTwoBodySphericalHarmonicTorque > torqueModel,
        const std::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
        const std::shared_ptr< simulation_setup::Body > bodyExertingTorque,
        const std::string& acceleratedBody,
        const std::string& acceleratingBody ):
    TorquePartial( acceleratedBody, acceleratingBody, basic_astrodynamics::full_two_body_spherical_harmonic_gravitational_torque ),
    torqueModel_( torqueModel ),
    bodyUndergoingTorqueObject_( bodyUndergoingTorque ),
    bodyExertingTorqueObject_( bodyExertingTorque ),
    currentPartialWrtQuaternionOfBodyUndergoingTorque_( Eigen::Matrix< double, 3, 4 >::Zero( ) ),
    currentPartialWrtQuaternionOfBodyExertingTorque_( Eigen::Matrix< double, 3, 4 >::Zero( ) ),
    currentPartialWrtPositionOfBodyUndergoingTorque_( Eigen::Matrix3d::Zero( ) ),
    currentPartialWrtPositionOfBodyExertingTorque_( Eigen::Matrix3d::Zero( ) ),
    quaternionPerturbation_( Eigen::Vector4d::Constant( 1.0E-9 ) ),
    positionPerturbation_( Eigen::Vector3d::Constant( 1.0 ) ),
    sphericalHarmonicCoefficientPerturbation_( 1.0E-8 )
{
    if( torqueModel_ == nullptr )
    {
        throw std::runtime_error(
                "Error when creating FullTwoBodySphericalHarmonicGravitationalTorquePartial, torque model is nullptr." );
    }
    if( bodyUndergoingTorqueObject_ == nullptr || bodyExertingTorqueObject_ == nullptr )
    {
        throw std::runtime_error(
                "Error when creating FullTwoBodySphericalHarmonicGravitationalTorquePartial, one or more body objects are nullptr." );
    }
}

std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
FullTwoBodySphericalHarmonicGravitationalTorquePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > )
{
    return std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );
}

std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
FullTwoBodySphericalHarmonicGravitationalTorquePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    using namespace estimatable_parameters;

    if( ( parameter->getParameterName( ).first == spherical_harmonics_cosine_coefficient_block ||
          parameter->getParameterName( ).first == spherical_harmonics_sine_coefficient_block ) &&
        ( parameter->getParameterName( ).second.first == bodyUndergoingTorque_ ||
          parameter->getParameterName( ).second.first == bodyExertingTorque_ ) )
    {
        return std::make_pair(
                std::bind( &FullTwoBodySphericalHarmonicGravitationalTorquePartial::wrtSphericalHarmonicCoefficientParameter,
                           this,
                           std::placeholders::_1,
                           parameter,
                           Eigen::VectorXd::Constant( parameter->getParameterSize( ), sphericalHarmonicCoefficientPerturbation_ ) ),
                parameter->getParameterSize( ) );
    }

    return std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );
}

void FullTwoBodySphericalHarmonicGravitationalTorquePartial::wrtSphericalHarmonicCoefficientParameter(
        Eigen::MatrixXd& partialMatrix,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > >& parameter,
        const Eigen::VectorXd& perturbation )
{
    partialMatrix = calculateTorqueWrtParameterPartials( parameter, torqueModel_, perturbation, emptyFunction, currentTime_ );
}

void FullTwoBodySphericalHarmonicGravitationalTorquePartial::wrtOrientationOfAcceleratedBody(
        Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const bool addContribution,
        const int startRow,
        const int startColumn )
{
    if( addContribution )
    {
        partialMatrix.block( startRow, startColumn, 3, 4 ) += currentPartialWrtQuaternionOfBodyUndergoingTorque_;
    }
    else
    {
        partialMatrix.block( startRow, startColumn, 3, 4 ) -= currentPartialWrtQuaternionOfBodyUndergoingTorque_;
    }
}

void FullTwoBodySphericalHarmonicGravitationalTorquePartial::wrtOrientationOfAcceleratingBody(
        Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const bool addContribution,
        const int startRow,
        const int startColumn )
{
    if( addContribution )
    {
        partialMatrix.block( startRow, startColumn, 3, 4 ) += currentPartialWrtQuaternionOfBodyExertingTorque_;
    }
    else
    {
        partialMatrix.block( startRow, startColumn, 3, 4 ) -= currentPartialWrtQuaternionOfBodyExertingTorque_;
    }
}

bool FullTwoBodySphericalHarmonicGravitationalTorquePartial::isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
        const std::pair< std::string, std::string >& stateReferencePoint,
        const propagators::IntegratedStateType integratedStateType )
{
    return ( integratedStateType == propagators::translational_state &&
             ( stateReferencePoint.first == bodyUndergoingTorque_ || stateReferencePoint.first == bodyExertingTorque_ ) );
}

void FullTwoBodySphericalHarmonicGravitationalTorquePartial::wrtNonRotationalStateOfAdditionalBody(
        Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const std::pair< std::string, std::string >& stateReferencePoint,
        const propagators::IntegratedStateType integratedStateType )
{
    if( integratedStateType == propagators::translational_state )
    {
        if( stateReferencePoint.first == bodyUndergoingTorque_ )
        {
            partialMatrix.block( 0, 0, 3, 3 ) += currentPartialWrtPositionOfBodyUndergoingTorque_;
        }
        else if( stateReferencePoint.first == bodyExertingTorque_ )
        {
            partialMatrix.block( 0, 0, 3, 3 ) += currentPartialWrtPositionOfBodyExertingTorque_;
        }
    }
}

void FullTwoBodySphericalHarmonicGravitationalTorquePartial::update( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        torqueModel_->updateMembers( currentTime );

        const Eigen::Vector7d rotationalStateOfBodyUndergoingTorque = bodyUndergoingTorqueObject_->getRotationalStateVector( );
        const Eigen::Vector7d rotationalStateOfBodyExertingTorque = bodyExertingTorqueObject_->getRotationalStateVector( );
        const Eigen::Vector6d translationalStateOfBodyUndergoingTorque = bodyUndergoingTorqueObject_->getState( );
        const Eigen::Vector6d translationalStateOfBodyExertingTorque = bodyExertingTorqueObject_->getState( );

        const std::function< void( Eigen::Vector7d ) > setRotationalStateOfBodyUndergoingTorque =
                std::bind( &simulation_setup::Body::setCurrentRotationalStateToLocalFrame,
                           bodyUndergoingTorqueObject_,
                           std::placeholders::_1 );
        const std::function< void( Eigen::Vector7d ) > setRotationalStateOfBodyExertingTorque =
                std::bind( &simulation_setup::Body::setCurrentRotationalStateToLocalFrame,
                           bodyExertingTorqueObject_,
                           std::placeholders::_1 );
        const std::function< void( Eigen::Vector6d ) > setTranslationalStateOfBodyUndergoingTorque =
                std::bind( &simulation_setup::Body::setState, bodyUndergoingTorqueObject_, std::placeholders::_1 );
        const std::function< void( Eigen::Vector6d ) > setTranslationalStateOfBodyExertingTorque =
                std::bind( &simulation_setup::Body::setState, bodyExertingTorqueObject_, std::placeholders::_1 );

        currentPartialWrtQuaternionOfBodyUndergoingTorque_ = calculateTorqueWrtRotationalStatePartials(
                setRotationalStateOfBodyUndergoingTorque,
                torqueModel_,
                rotationalStateOfBodyUndergoingTorque,
                quaternionPerturbation_,
                0,
                4,
                emptyFunction,
                currentTime );

        currentPartialWrtQuaternionOfBodyExertingTorque_ = calculateTorqueWrtRotationalStatePartials(
                setRotationalStateOfBodyExertingTorque,
                torqueModel_,
                rotationalStateOfBodyExertingTorque,
                quaternionPerturbation_,
                0,
                4,
                emptyFunction,
                currentTime );

        currentPartialWrtPositionOfBodyUndergoingTorque_ = calculateTorqueWrtTranslationalStatePartials(
                setTranslationalStateOfBodyUndergoingTorque,
                torqueModel_,
                translationalStateOfBodyUndergoingTorque,
                positionPerturbation_,
                0,
                emptyFunction,
                currentTime );

        currentPartialWrtPositionOfBodyExertingTorque_ = calculateTorqueWrtTranslationalStatePartials(
                setTranslationalStateOfBodyExertingTorque,
                torqueModel_,
                translationalStateOfBodyExertingTorque,
                positionPerturbation_,
                0,
                emptyFunction,
                currentTime );

        currentTime_ = currentTime;
    }
}

}  // namespace acceleration_partials

}  // namespace tudat

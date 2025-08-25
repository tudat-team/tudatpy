/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/acceleration_partials/aerodynamicAccelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

void AerodynamicAccelerationPartial::computeAerodynamicAccelerationWrtDragComponent(
    Eigen::MatrixXd& partial )
{
    Eigen::Quaterniond rotationToInertialFrame =
                flightConditions_->getAerodynamicAngleCalculator( )->getRotationQuaternionBetweenFrames(
                        reference_frames::aerodynamic_frame, reference_frames::inertial_frame );

    Eigen::Vector3d currentDragComponentPartial = Eigen::Vector3d::Zero( );
    Eigen::Vector3d unscaledAcceleration = aerodynamicAcceleration_->getCurrentUnscaledAccelerationInAerodynamicFrame( );
    currentDragComponentPartial( 0 ) = unscaledAcceleration( 0 );
    partial = rotationToInertialFrame * currentDragComponentPartial;
};

void AerodynamicAccelerationPartial::computeAerodynamicAccelerationWrtSideComponent(
    Eigen::MatrixXd& partial )
{
    Eigen::Quaterniond rotationToInertialFrame =
                flightConditions_->getAerodynamicAngleCalculator( )->getRotationQuaternionBetweenFrames(
                        reference_frames::aerodynamic_frame, reference_frames::inertial_frame );

    Eigen::Vector3d currentSideComponentPartial = Eigen::Vector3d::Zero( );
    Eigen::Vector3d unscaledAcceleration = aerodynamicAcceleration_->getCurrentUnscaledAccelerationInAerodynamicFrame( );
    currentSideComponentPartial( 1 ) = unscaledAcceleration( 1 );
    partial = rotationToInertialFrame * currentSideComponentPartial;
};

void AerodynamicAccelerationPartial::computeAerodynamicAccelerationWrtLiftComponent(
    Eigen::MatrixXd& partial )
{
    Eigen::Quaterniond rotationToInertialFrame =
                flightConditions_->getAerodynamicAngleCalculator( )->getRotationQuaternionBetweenFrames(
                        reference_frames::aerodynamic_frame, reference_frames::inertial_frame );

    Eigen::Vector3d currentLiftComponentPartial = Eigen::Vector3d::Zero( );
    Eigen::Vector3d unscaledAcceleration = aerodynamicAcceleration_->getCurrentUnscaledAccelerationInAerodynamicFrame( );
    currentLiftComponentPartial( 2 ) = unscaledAcceleration( 2 );
    partial = rotationToInertialFrame * currentLiftComponentPartial;  
};

//! Function for updating partial w.r.t. the bodies' positions
void AerodynamicAccelerationPartial::update( const double currentTime )
{
    Eigen::Vector6d nominalState = vehicleStateGetFunction_( );
    Eigen::Vector6d perturbedState;

    // Compute state partial by numerical difference
    Eigen::Vector3d upperturbedAcceleration, downperturbedAcceleration;
    for( unsigned int i = 0; i < 6; i++ )
    {
        // Perturb state upwards
        perturbedState = nominalState;
        perturbedState( i ) += bodyStatePerturbations_( i );

        // Update environment/acceleration to perturbed state.
        flightConditions_->resetCurrentTime( );
        aerodynamicAcceleration_->resetCurrentTime( );
        vehicleStateSetFunction_( perturbedState );
        flightConditions_->updateConditions( currentTime );
        aerodynamicAcceleration_->updateMembers( currentTime );

        // Retrieve perturbed acceleration.
        upperturbedAcceleration = aerodynamicAcceleration_->getAcceleration( );

        // Perturb state downwards
        perturbedState = nominalState;
        perturbedState( i ) -= bodyStatePerturbations_( i );

        // Update environment/acceleration to perturbed state.
        flightConditions_->resetCurrentTime( );
        aerodynamicAcceleration_->resetCurrentTime( );
        vehicleStateSetFunction_( perturbedState );
        flightConditions_->updateConditions( currentTime );
        aerodynamicAcceleration_->updateMembers( currentTime );

        // Retrieve perturbed acceleration.
        downperturbedAcceleration = aerodynamicAcceleration_->getAcceleration( );

        // Compute partial
        currentAccelerationStatePartials_.block( 0, i, 3, 1 ) =
                ( upperturbedAcceleration - downperturbedAcceleration ) / ( 2.0 * bodyStatePerturbations_( i ) );
    }

    // Reset environment/acceleration mode to nominal conditions
    flightConditions_->resetCurrentTime( );
    aerodynamicAcceleration_->resetCurrentTime( );

    vehicleStateSetFunction_( nominalState );
    flightConditions_->updateConditions( currentTime );
    aerodynamicAcceleration_->updateMembers( currentTime );

    currentTime_ = currentTime;
}

std::pair< std::function< void( Eigen::MatrixXd& ) >, int > AerodynamicAccelerationPartial::getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfColumns = 0;
    if( parameter->getParameterName( ).second.first == acceleratedBody_ )
    {
        switch( parameter->getParameterName( ).first )
        {
            case estimatable_parameters::constant_drag_coefficient:
            {
                partialFunction = std::bind(
                    &AerodynamicAccelerationPartial::computeAccelerationPartialWrtCurrentDragCoefficient, this, std::placeholders::_1 );
                numberOfColumns = 1;
                break;
            }
            case estimatable_parameters::drag_component_scaling_factor:
            {
                partialFunction = std::bind(
                    &AerodynamicAccelerationPartial::computeAerodynamicAccelerationWrtDragComponent, 
                    this, std::placeholders::_1 );
                numberOfColumns = 1;
                break;
            }
            case estimatable_parameters::side_component_scaling_factor:
            {
                partialFunction = std::bind(
                    &AerodynamicAccelerationPartial::computeAerodynamicAccelerationWrtSideComponent, 
                    this, std::placeholders::_1 );
                numberOfColumns = 1;
                break;
            }
            case estimatable_parameters::lift_component_scaling_factor:
            {
                partialFunction = std::bind(
                    &AerodynamicAccelerationPartial::computeAerodynamicAccelerationWrtLiftComponent, 
                    this, std::placeholders::_1 );
                numberOfColumns = 1;
                break;
            }
            default:
                break;
        }
    }

    if( numberOfColumns == 0 )
    {
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > basePartialFunctionPair =
                this->getParameterPartialFunctionAccelerationBase( parameter );
        if( basePartialFunctionPair.second != 0 )
        {
            partialFunction = basePartialFunctionPair.first;
            numberOfColumns = basePartialFunctionPair.second;
        }
    }

    return std::make_pair( partialFunction, numberOfColumns );
}

}  // namespace acceleration_partials

}  // namespace tudat

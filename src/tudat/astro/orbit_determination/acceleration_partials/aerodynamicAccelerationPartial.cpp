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

//! Function to compute the partial derivative of the acceleration w.r.t. the current drag coefficient
/*!
 * Function to compute the partial derivative of the acceleration w.r.t. the current drag coefficient
 * \param accelerationPartial Derivative of acceleration w.r.t. the current drag coefficient (returned by reference).
 */

void AerodynamicAccelerationPartial::computeAccelerationPartialWrtCurrentDragCoefficient( Eigen::MatrixXd& accelerationPartial )
{
    Eigen::Quaterniond rotationToInertialFrame = flightConditions_->getAerodynamicAngleCalculator( )->getRotationQuaternionBetweenFrames(
            reference_frames::aerodynamic_frame, reference_frames::inertial_frame );

    double currentAirspeed = flightConditions_->getCurrentAirspeed( );
    accelerationPartial = rotationToInertialFrame * Eigen::Vector3d::UnitX( ) *
            ( -0.5 * flightConditions_->getCurrentDensity( ) * currentAirspeed * currentAirspeed *
              flightConditions_->getAerodynamicCoefficientInterface( )->getReferenceArea( ) ) /
            aerodynamicAcceleration_->getCurrentMass( );
}

//! Function to compute the partial derivative of the acceleration w.r.t. the arc-wise constant drag coefficient
/*!
 * Function to compute the partial derivative of the acceleration w.r.t. the arc-wise constant drag coefficient
 * \param accelerationPartial Derivative of acceleration w.r.t. arc-wise constant drag coefficient (returned by reference).
 * \param parameter Parameter object containing information on arcwise drag coefficient that is to be estimated
 */

void AerodynamicAccelerationPartial::computeAccelerationPartialWrtArcwiseDragCoefficient(
        Eigen::MatrixXd& accelerationPartial,
        const std::shared_ptr< estimatable_parameters::ArcWiseConstantDragCoefficient > parameter )
{
    // Get partial w.r.t. rdrag coefficient
    Eigen::MatrixXd partialWrtSingleParameter = Eigen::Vector3d::Zero( );
    this->computeAccelerationPartialWrtCurrentDragCoefficient( partialWrtSingleParameter );

    // Retrieve current arc
    std::shared_ptr< interpolators::LookUpScheme< double > > currentArcIndexLookUp = parameter->getArcTimeLookupScheme( );
    accelerationPartial.setZero( 3, parameter->getNumberOfArcs( ) );
    if( currentArcIndexLookUp->getMinimumValue( ) <= currentTime_ )
    {
        int currentArc = currentArcIndexLookUp->findNearestLowerNeighbour( currentTime_ );

        if( currentArc >= accelerationPartial.cols( ) )
        {
            throw std::runtime_error( "Error when getting arc-wise radiation pressure coefficient partials, data not consistent" );
        }

        // Set partial
        accelerationPartial.block( 0, currentArc, 3, 1 ) = partialWrtSingleParameter;
    }
}

void AerodynamicAccelerationPartial::computeAccelerationPartialWrtDragComponent( Eigen::MatrixXd& partial )
{
    Eigen::Quaterniond rotationToInertialFrame = flightConditions_->getAerodynamicAngleCalculator( )->getRotationQuaternionBetweenFrames(
            reference_frames::aerodynamic_frame, reference_frames::inertial_frame );

    Eigen::Vector3d currentDragComponentPartial = Eigen::Vector3d::Zero( );
    Eigen::Vector3d unscaledAcceleration = aerodynamicAcceleration_->getCurrentUnscaledAccelerationInAerodynamicFrame( );
    currentDragComponentPartial( 0 ) = unscaledAcceleration( 0 );
    partial = rotationToInertialFrame * currentDragComponentPartial;
};

void AerodynamicAccelerationPartial::computeAccelerationPartialWrtSideComponent( Eigen::MatrixXd& partial )
{
    Eigen::Quaterniond rotationToInertialFrame = flightConditions_->getAerodynamicAngleCalculator( )->getRotationQuaternionBetweenFrames(
            reference_frames::aerodynamic_frame, reference_frames::inertial_frame );

    Eigen::Vector3d currentSideComponentPartial = Eigen::Vector3d::Zero( );
    Eigen::Vector3d unscaledAcceleration = aerodynamicAcceleration_->getCurrentUnscaledAccelerationInAerodynamicFrame( );
    currentSideComponentPartial( 1 ) = unscaledAcceleration( 1 );
    partial = rotationToInertialFrame * currentSideComponentPartial;
};

void AerodynamicAccelerationPartial::computeAccelerationPartialWrtLiftComponent( Eigen::MatrixXd& partial )
{
    Eigen::Quaterniond rotationToInertialFrame = flightConditions_->getAerodynamicAngleCalculator( )->getRotationQuaternionBetweenFrames(
            reference_frames::aerodynamic_frame, reference_frames::inertial_frame );

    Eigen::Vector3d currentLiftComponentPartial = Eigen::Vector3d::Zero( );
    Eigen::Vector3d unscaledAcceleration = aerodynamicAcceleration_->getCurrentUnscaledAccelerationInAerodynamicFrame( );
    currentLiftComponentPartial( 2 ) = unscaledAcceleration( 2 );
    partial = rotationToInertialFrame * currentLiftComponentPartial;
};

//! Function to compute the partial derivative of the acceleration w.r.t. an aerodynamic component scaling factor
/*!
 * Function to compute the partial derivative of the acceleration w.r.t. an aerodynamic component scaling factor
 * \param accelerationPartial Derivative of acceleration w.r.t. an aerodynamic component scaling factor (returned by reference).
 */
void AerodynamicAccelerationPartial::computeAccelerationPartialWrtArcWiseAerodynamicScalingCofficient(
        Eigen::MatrixXd& accelerationPartial,
        const std::shared_ptr< estimatable_parameters::ArcWiseAerodynamicScalingFactor > parameter )
{
    // Get partial w.r.t. aerodynamic component coefficient
    Eigen::MatrixXd partialWrtSingleParameter = Eigen::Vector3d::Zero( );
    switch( parameter->getParameterName( ).first )
    {
        case estimatable_parameters::arc_wise_drag_component_scaling_factor: {
            this->computeAccelerationPartialWrtDragComponent( partialWrtSingleParameter );
            break;
        }
        case estimatable_parameters::arc_wise_side_component_scaling_factor: {
            this->computeAccelerationPartialWrtSideComponent( partialWrtSingleParameter );
            break;
        }
        case estimatable_parameters::arc_wise_lift_component_scaling_factor: {
            this->computeAccelerationPartialWrtLiftComponent( partialWrtSingleParameter );
            break;
        }

        default:
            break;
    }

    // Retrieve current arc
    std::shared_ptr< interpolators::LookUpScheme< double > > currentArcIndexLookUp = parameter->getArcTimeLookupScheme( );
    accelerationPartial.setZero( 3, parameter->getNumberOfArcs( ) );
    if( currentArcIndexLookUp->getMinimumValue( ) <= currentTime_ )
    {
        int currentArc = currentArcIndexLookUp->findNearestLowerNeighbour( currentTime_ );

        if( currentArc >= accelerationPartial.cols( ) )
        {
            throw std::runtime_error( "Error when getting arc-wise aerodynamic component scaling partials, data not consistent" );
        }

        // Set partial
        accelerationPartial.block( 0, currentArc, 3, 1 ) = partialWrtSingleParameter;
    }
}



//! Function to compute the partial derivative of the acceleration w.r.t. the current atmospheric density
/*!
 * Function to compute the partial derivative of the acceleration w.r.t. the current atmospheric density
 * \param accelerationPartial Derivative of acceleration by reference.
 */
void AerodynamicAccelerationPartial::computeAccelerationPartialWrtCurrentDensity( Eigen::MatrixXd& accelerationPartial )
{

    Eigen::Quaterniond rotationToInertialFrame = flightConditions_->getAerodynamicAngleCalculator( )->getRotationQuaternionBetweenFrames(
        reference_frames::aerodynamic_frame, reference_frames::inertial_frame );

    // retrieve from flight conditions
    Eigen::Vector3d currentForceCoefficients = flightConditions_->getAerodynamicCoefficientInterface( )->getCurrentForceCoefficients(  );
    double const currentAirspeed = flightConditions_->getCurrentAirspeed( );
    double const referenceArea = flightConditions_->getAerodynamicCoefficientInterface( )->getReferenceArea( );
    double const currentMass = aerodynamicAcceleration_->getCurrentMass( );

    // 1/2 . Ci . V2 . A / M
    accelerationPartial = rotationToInertialFrame * (0.5 * currentAirspeed * currentAirspeed * referenceArea / currentMass * currentForceCoefficients);

}



//! Function to compute the partial derivative of the acceleration w.r.t. base density
/*!
 * Function to compute the partial derivative of the acceleration w.r.t. base density
 * \param accelerationPartial Derivative of acceleration w.r.t. base density (by reference)
 */

void AerodynamicAccelerationPartial::computeAccelerationPartialWrtExponentialAtmosphereBaseDensity( Eigen::MatrixXd& accelerationPartial )
{
    // dA/dRho:
    this->computeAccelerationPartialWrtCurrentDensity(accelerationPartial);
    // dRho/dRho_0
    double partialCurrentDensityWrtBaseDensity = std::exp( -flightConditions_->getCurrentAltitude(  ) / exponentialAtmosphere_->getScaleHeight( ) );

    accelerationPartial *= partialCurrentDensityWrtBaseDensity;
}


//! Function to compute the partial derivative of the acceleration w.r.t. scale height
/*!
 * Function to compute the partial derivative of the acceleration w.r.t. scale height
 * \param accelerationPartial Derivative of acceleration w.r.t. scale height (by reference)
 */

void AerodynamicAccelerationPartial::computeAccelerationPartialWrtExponentialAtmosphereScaleHeight( Eigen::MatrixXd& accelerationPartial )
{
    // dA/dRho:
    this->computeAccelerationPartialWrtCurrentDensity(accelerationPartial);
    // dRho/dH
    double scaleHeight = exponentialAtmosphere_->getScaleHeight( );
    double expTerm = std::exp( -flightConditions_->getCurrentAltitude(  ) / scaleHeight );
    double baseTerm = exponentialAtmosphere_->getBaseDensity( ) * flightConditions_->getCurrentAltitude(  ) / scaleHeight / scaleHeight;
    double partialCurrentDensityWrtScaleHeight = baseTerm*expTerm;

    accelerationPartial *= partialCurrentDensityWrtScaleHeight;
}


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

// overwrite base class (essential!) function, which evaluates all acceleration model parameter (double) partials:
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > AerodynamicAccelerationPartial::getParameterPartialFunctionDerivedAcceleration(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfColumns = 0;
    if( parameter->getParameterName( ).second.first == acceleratedBody_ )
    {
        switch( parameter->getParameterName( ).first )
        {
            case estimatable_parameters::constant_drag_coefficient: {
                partialFunction = std::bind(
                        &AerodynamicAccelerationPartial::computeAccelerationPartialWrtCurrentDragCoefficient, this, std::placeholders::_1 );
                numberOfColumns = 1;
                break;
            }
            case estimatable_parameters::drag_component_scaling_factor: {
                partialFunction = std::bind(
                        &AerodynamicAccelerationPartial::computeAccelerationPartialWrtDragComponent, this, std::placeholders::_1 );
                numberOfColumns = 1;
                break;
            }
            case estimatable_parameters::side_component_scaling_factor: {
                partialFunction = std::bind(
                        &AerodynamicAccelerationPartial::computeAccelerationPartialWrtSideComponent, this, std::placeholders::_1 );
                numberOfColumns = 1;
                break;
            }
            case estimatable_parameters::lift_component_scaling_factor: {
                partialFunction = std::bind(
                        &AerodynamicAccelerationPartial::computeAccelerationPartialWrtLiftComponent, this, std::placeholders::_1 );
                numberOfColumns = 1;
                break;
            }
            default:
                break;
        }
    }

    return std::make_pair( partialFunction, numberOfColumns );
}

// overwrite base class (essential!) function, which evaluates all acceleration model parameter (VectorXd) partials:
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > AerodynamicAccelerationPartial::getParameterPartialFunctionDerivedAcceleration(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )

{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfColumns = 0;
    if( parameter->getParameterName( ).second.first == acceleratedBody_ )
    {
        switch( parameter->getParameterName( ).first )
        {
            case estimatable_parameters::arc_wise_constant_drag_coefficient: {
                if( std::dynamic_pointer_cast< estimatable_parameters::ArcWiseConstantDragCoefficient >( parameter ) != nullptr )
                {
                    partialFunction =
                            std::bind( &AerodynamicAccelerationPartial::computeAccelerationPartialWrtArcwiseDragCoefficient,
                                       this,
                                       std::placeholders::_1,
                                       std::dynamic_pointer_cast< estimatable_parameters::ArcWiseConstantDragCoefficient >( parameter ) );
                    numberOfColumns = parameter->getParameterSize( );
                }
                else
                {
                    throw std::runtime_error( "Error when making radiation drag, arcwise drag parameter not consistent" );
                }
                break;
            }
            case estimatable_parameters::arc_wise_drag_component_scaling_factor:
            case estimatable_parameters::arc_wise_side_component_scaling_factor:
            case estimatable_parameters::arc_wise_lift_component_scaling_factor: {
                if( std::dynamic_pointer_cast< estimatable_parameters::ArcWiseAerodynamicScalingFactor >( parameter ) != nullptr )
                {
                    partialFunction =
                            std::bind( &AerodynamicAccelerationPartial::computeAccelerationPartialWrtArcWiseAerodynamicScalingCofficient,
                                       this,
                                       std::placeholders::_1,
                                       std::dynamic_pointer_cast< estimatable_parameters::ArcWiseAerodynamicScalingFactor >( parameter ) );
                    numberOfColumns = parameter->getParameterSize( );
                }
                else
                {
                    throw std::runtime_error( "Error when making radiation drag, arcwise drag parameter not consistent" );
                }
                break;
            }
            default:
                break;
        }
    }

    return std::make_pair( partialFunction, numberOfColumns );
}

}  // namespace acceleration_partials

}  // namespace tudat

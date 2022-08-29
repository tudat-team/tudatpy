/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/acceleration_partials/thrustAccelerationPartial.h"


namespace tudat
{

namespace acceleration_partials
{

ThrustAccelerationPartial::ThrustAccelerationPartial(
        const std::shared_ptr< propulsion::ThrustAcceleration > thrustAcceleration,
        const std::string acceleratedBody,
        const std::map< std::pair< estimatable_parameters::EstimatebleParametersEnum, std::string >,
                std::shared_ptr< observation_partials::RotationMatrixPartial > >& rotationMatrixPartials  ):
    AccelerationPartial( acceleratedBody, acceleratedBody, basic_astrodynamics::thrust_acceleration ),
    thrustAcceleration_( thrustAcceleration ),
    rotationMatrixPartials_( rotationMatrixPartials ),
    isAccelerationDependentOnTranslationalState_( false )
{
    thrustAcceleration_->setSaveThrustContributions( true );
    thrustSources_ = thrustAcceleration_-> getThrustSources( );

    isAccelerationDependentOnMass_ = false;
    std::vector< std::shared_ptr< system_models::EngineModel > > thrustSources =
            thrustAcceleration_->getThrustSources( );

    for( unsigned int i = 0; i < thrustSources.size( ); i++ )
    {
        if( thrustSources.at( i )->getThrustMagnitudeWrapper( )->modelIsForceBased( ) )
        {
            isAccelerationDependentOnMass_ = true;
            massDependentThrustSources_.push_back( i );
        }

        if( std::dynamic_pointer_cast< propulsion::ParameterizedThrustMagnitudeWrapper >(
                    thrustSources.at( i )->getThrustMagnitudeWrapper( ) ) )
        {
            std::cerr<<"Warning, engine "<<thrustSources.at( i )->getEngineName( )<<
                       " has a ParameterizedThrustMagnitudeWrapper, thrust acceleration partials will no be calculated properly"<<std::endl;
        }

    }


    if( std::dynamic_pointer_cast< propulsion::DirectThrustDirectionCalculator >(
                thrustAcceleration->getThrustDirectionCalculator( ) ) != nullptr )
    {
        isRotationDirectionBased_ = true;
    }
    else
    {
        isRotationDirectionBased_ = false;
    }

}


bool ThrustAccelerationPartial::isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
        const std::pair< std::string, std::string >& stateReferencePoint,
        const propagators::IntegratedStateType integratedStateType )
{
    bool isStateDependent = false;
    if( integratedStateType == propagators::body_mass_state && stateReferencePoint.first == acceleratedBody_ )
    {
        isStateDependent = isAccelerationDependentOnMass_;
    }

    if( integratedStateType == propagators::rotational_state && stateReferencePoint.first == acceleratedBody_ )
    {
        isStateDependent = true;
    }
    return isStateDependent;
}

//! Function to calculate an acceleration partial wrt a rotational parameter.
void ThrustAccelerationPartial::wrtRotationModelParameter(
        Eigen::MatrixXd& accelerationPartial,
        const estimatable_parameters::EstimatebleParametersEnum parameterType,
        const std::string& secondaryIdentifier )
{

    if( rotationMatrixPartials_.count( std::make_pair( parameterType, secondaryIdentifier ) ) == 0 )
    {
        throw std::runtime_error( "Error when calculating thrust parial w.r.t. rotation matrix paramater, calculator object not found for " +
                                  std::to_string( parameterType ) + ", " + secondaryIdentifier );
    }

    // Get rotation matrix partial(s) wrt requested parameter
    std::vector< Eigen::Matrix3d > rotationMatrixPartials =
            rotationMatrixPartials_.at( std::make_pair( parameterType, secondaryIdentifier ) )->
            calculatePartialOfRotationMatrixToBaseFrameWrParameter( currentTime_ );

    // Compute total body-fixed thrust
    Eigen::Vector3d currentBodyFixedThrust = Eigen::Vector3d::Zero( );
    for( unsigned int j = 0; j < thrustAcceleration_->getThrustSources( ).size( ); j++ )
    {
        currentBodyFixedThrust += thrustSources_.at( j )->getBodyFixedThrustDirection( ) *
                thrustSources_.at( j )->getCurrentThrustAcceleration( thrustAcceleration_->getCurrentBodyMass( ) );

    }

    // Compute derivate w.r.t. each parameter entry
    for( unsigned int i = 0; i < rotationMatrixPartials.size( ); i++ )
    {
        accelerationPartial.block( 0, i, 3, 1 ) += rotationMatrixPartials[ i ] * currentBodyFixedThrust;

    }
}

void ThrustAccelerationPartial::wrtBodyMass( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                             const bool addContribution )
{
    partialMatrix.setZero( );
    for( unsigned int i = 0; i < massDependentThrustSources_.size( ); i++ )
    {
        partialMatrix.block( 0, 0, 3, 1 ) +=
                ( addContribution ? 1.0 : -1.0 ) * thrustAcceleration_->getCurrentThrustAccelerationContribution(
                    massDependentThrustSources_.at( i ) );

    }
    partialMatrix.block( 0, 0, 3, 1 ) /= -thrustAcceleration_->getCurrentBodyMass( );
}

void ThrustAccelerationPartial::wrtThrustMagnitude(
        Eigen::MatrixXd& partialMatrix,
        const int engineIndex )
{
    partialMatrix.setZero( );
    partialMatrix.block( 0, 0, 3, 1 ) += thrustAcceleration_->getCurrentThrustAccelerationContribution(
                engineIndex ) / thrustSources_.at( engineIndex )->getCurrentThrust( thrustAcceleration_->getCurrentBodyMass( ) );
}

void ThrustAccelerationPartial::wrtNonTranslationalStateOfAdditionalBody(
        Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const std::pair< std::string, std::string >& stateReferencePoint,
        const propagators::IntegratedStateType integratedStateType,
        const bool addContribution )
{
    // If partial is w.r.t. body mass, iterate over all mass-dependent sources, and compute contributions
    if( integratedStateType == propagators::body_mass_state && stateReferencePoint.first == acceleratedBody_ )
    {
        wrtBodyMass( partialMatrix, addContribution );
    }

    // If partial is w.r.t. rotational state, call corresponding rotation matrix partial
    if( integratedStateType == propagators::rotational_state && stateReferencePoint.first == acceleratedBody_ )
    {
        Eigen::MatrixXd tempMatrix = Eigen::MatrixXd::Zero( 3, 7 );
        wrtRotationModelParameter( tempMatrix, estimatable_parameters::initial_rotational_body_state, "" );
        partialMatrix.block( 0, 0, 3, 7 ) = ( addContribution ? 1.0 : -1.0 ) * tempMatrix;
    }
}

//! Constructor
MomentumWheelDesaturationPartial::MomentumWheelDesaturationPartial(
        const std::shared_ptr< propulsion::MomentumWheelDesaturationThrustAcceleration > thrustAcceleration,
        const std::string acceleratedBody ):
    AccelerationPartial( acceleratedBody, acceleratedBody, basic_astrodynamics::momentum_wheel_desaturation_acceleration ),
    thrustAcceleration_( thrustAcceleration ){ }

//! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
MomentumWheelDesaturationPartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )

{
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionPair;

    // Check dependencies.
    if( parameter->getParameterName( ).first == estimatable_parameters::desaturation_delta_v_values )
    {
        // If parameter is desaturation deltaV values, check and create dependency function .
        partialFunctionPair = std::make_pair(
                    std::bind( &MomentumWheelDesaturationPartial::wrtDesaturationDeltaVValues, this, std::placeholders::_1 ),
                    parameter->getParameterSize( ) );
    }
    else
    {
        partialFunctionPair = std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );
    }

    return partialFunctionPair;
}


//! Function to compute the partial derivative w.r.t. the deltaV values of the momentum desaturation maneuvers
void MomentumWheelDesaturationPartial::wrtDesaturationDeltaVValues( Eigen::MatrixXd& accelerationPartial )
{
    // Compute partials.
    accelerationPartial.setZero( );
    accelerationPartial.block( 0, 3 * thrustAcceleration_->getCurrentNearestTimeIndex( ), 3, 3 ) =
          Eigen::Matrix3d::Identity( ) * thrustAcceleration_->getCurrentThrustMultiplier( ) /
              ( thrustAcceleration_->getTotalManeuverTime( ) - thrustAcceleration_->getManeuverRiseTime( ) );
}

}

}

/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/acceleration_partials/fullRadiationPressureAccelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

RadiationPressureAccelerationPartial::RadiationPressureAccelerationPartial(
    const std::string acceleratedBody,
    const std::string acceleratingBody,
    const std::shared_ptr< electromagnetism::PaneledSourceRadiationPressureAcceleration > radiationPressureAcceleration,
    const std::shared_ptr< estimatable_parameters::CustomSingleAccelerationPartialCalculatorSet > customAccelerationPartialSet ):
    AccelerationPartial( acceleratedBody, acceleratingBody, basic_astrodynamics::custom_acceleration ),
    radiationPressureAcceleration_( radiationPressureAcceleration ),
    customAccelerationPartialSet_( customAccelerationPartialSet )
{
    currentPartialWrtUndergoingState_.setZero( );
    currentPartialWrtExertingState_.setZero( );

    estimatable_parameters::EstimatebleParameterIdentifier undergoingBodyIdentifier =
        std::make_pair( estimatable_parameters::initial_body_state, std::make_pair( acceleratedBody, "" ) );
    if( customAccelerationPartialSet->customInitialStatePartials_.count( undergoingBodyIdentifier ) > 0 )
    {
        bodyUndergoingPositionPartial_ = customAccelerationPartialSet->customInitialStatePartials_.at( undergoingBodyIdentifier );
    }

    estimatable_parameters::EstimatebleParameterIdentifier exertingBodyIdentifier =
        std::make_pair( estimatable_parameters::initial_body_state, std::make_pair( acceleratingBody, "" ) );
    if( customAccelerationPartialSet->customInitialStatePartials_.count( exertingBodyIdentifier ) > 0 && ( exertingBodyIdentifier != undergoingBodyIdentifier ) )
    {
        bodyExertingPositionPartial_ = customAccelerationPartialSet->customInitialStatePartials_.at( exertingBodyIdentifier );
    }
}

void RadiationPressureAccelerationPartial::update( const double currentTime )
{
    radiationPressureAcceleration_->updateMembers( currentTime );

    if( !( currentTime_ == currentTime ) )
    {
        if( bodyUndergoingPositionPartial_ != nullptr )
        {
            currentPartialWrtUndergoingState_ = bodyUndergoingPositionPartial_->computePartial(
                currentTime, radiationPressureAcceleration_->getAcceleration( ), radiationPressureAcceleration_ );
        }

        if( bodyExertingPositionPartial_ != nullptr )
        {
            currentPartialWrtExertingState_ = bodyExertingPositionPartial_->computePartial(
                currentTime, radiationPressureAcceleration_->getAcceleration( ), radiationPressureAcceleration_ );
        }

        currentTime_ = currentTime;
    }
}

void RadiationPressureAccelerationPartial::wrtRadiationPressureCoefficient(
    Eigen::MatrixXd& partial, std::shared_ptr< electromagnetism::CannonballRadiationPressureTargetModel > targetModel )
{
    if( targetModel->getCoefficient( ) == 0.0 )
    {
        throw std::runtime_error( "Error in full radiation pressure partial w.r.t. Cr, partial is only implemented for non-zero coefficient" );
    }
    partial = radiationPressureAcceleration_->getAcceleration( ) / targetModel->getCoefficient( );

}

}

}

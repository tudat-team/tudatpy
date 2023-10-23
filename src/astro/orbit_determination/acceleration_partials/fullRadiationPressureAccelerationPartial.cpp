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

    std::cout<<"Initial state partials "<<std::endl;
    for( auto it : customAccelerationPartialSet->customInitialStatePartials_ )
    {
        std::cout<<"Entries: "<<it.first.first<<" "<<it.first.second.first<<" "<<it.first.second.second<<std::endl;
    }
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

    std::cout<<currentTime<<" "<<currentTime_<<std::endl;
    if( !( currentTime_ == currentTime ) )
    {
        std::cout<<"Partials "<<bodyUndergoingPositionPartial_<<" "<<bodyExertingPositionPartial_<<std::endl;
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

}

}

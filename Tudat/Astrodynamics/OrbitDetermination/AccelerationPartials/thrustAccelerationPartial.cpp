/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/thrustAccelerationPartial.h"


namespace tudat
{

namespace acceleration_partials
{


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

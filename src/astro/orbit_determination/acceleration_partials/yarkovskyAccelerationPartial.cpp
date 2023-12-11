/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/acceleration_partials/yarkovskyAccelerationPartial.h"


namespace tudat
{

namespace acceleration_partials
{


Eigen::Matrix3d calculatePartialOfYarkovskyAccelerationWrtPositionOfAcceleratedBody(
    const Eigen::Vector6d& relativeState,
    const double yarkovskyParameter )
{
    const double currentDistance = relativeState.segment( 0, 3 ).norm( );
    const Eigen::Vector3d positionUnitVector = relativeState.segment( 0, 3 ).normalized( );
    const Eigen::Vector3d velocityUnitVector = relativeState.segment( 3, 3 ).normalized( );

    return -2.0 * yarkovskyParameter * physical_constants::ASTRONOMICAL_UNIT * physical_constants::ASTRONOMICAL_UNIT *
        ( velocityUnitVector * positionUnitVector.transpose( ) ) / ( currentDistance * currentDistance * currentDistance );
}

Eigen::Matrix3d calculatePartialOfYarkovskyAccelerationWrtVelocityOfAcceleratedBody(
    const Eigen::Vector6d& relativeState,
    const double yarkovskyParameter )
{
    const double currentDistance = relativeState.segment( 0, 3 ).norm( );
    const double currentSpeed = relativeState.segment( 3, 3 ).norm( );

    const Eigen::Vector3d velocityUnitVector = relativeState.segment( 3, 3 ).normalized( );

    return yarkovskyParameter * ( physical_constants::ASTRONOMICAL_UNIT * physical_constants::ASTRONOMICAL_UNIT ) / (currentDistance * currentDistance ) *
        ( Eigen::Matrix3d::Identity( ) - velocityUnitVector * velocityUnitVector.transpose( ) ) / currentSpeed;
}


//! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
YarkovskyAccelerationPartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )

{
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionPair;

    // Check dependencies.
    if( parameter->getParameterName( ).first == estimatable_parameters::yarkovsky_parameter &&
        parameter->getParameterName( ).second.first == this->getAcceleratedBody( ) )
    {
        // If parameter is gravitational parameter, check and create dependency function .
        partialFunctionPair =
            std::make_pair( std::bind( &YarkovskyAccelerationPartial::wrtYarkovskyParameter,
                       this, std::placeholders::_1 ), 1 );
    }
    else
    {
        partialFunctionPair = std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );
    }

    return partialFunctionPair;
}

//! Function to calculate central gravity partial w.r.t. central body gravitational parameter
void YarkovskyAccelerationPartial::wrtYarkovskyParameter( Eigen::MatrixXd& yarkovskyPartial )
{
    if( yarkovskyAcceleration_->getYarkovskyParameter( ) == 0.0 )
    {
        yarkovskyPartial = electromagnetism::computeYarkovskyAcceleration( 1.0, yarkovskyAcceleration_->getCurrentState( ) );
    }
    else
    {
        yarkovskyPartial = yarkovskyAcceleration_->getAcceleration( ) / yarkovskyAcceleration_->getYarkovskyParameter( );
    }
}



}

}

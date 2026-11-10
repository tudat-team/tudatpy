/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GROUNDSTATIONLINEARVELOCITY_H
#define TUDAT_GROUNDSTATIONLINEARVELOCITY_H

#include <memory>
#include <string>

#include "tudat/astro/ground_stations/groundStationState.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Body-fixed velocity (in m/s) of the linear motion model of a ground station.
class GroundStationLinearVelocity : public EstimatableParameter< Eigen::VectorXd >
{
public:
    GroundStationLinearVelocity( const std::shared_ptr< ground_stations::LinearStationMotionModel > motionModel,
                                 const std::string& associatedBody,
                                 const std::string& associatedStation ):
        EstimatableParameter< Eigen::VectorXd >( ground_station_linear_velocity, associatedBody, associatedStation ),
        motionModel_( motionModel )
    {}

    Eigen::VectorXd getParameterValue( )
    {
        return motionModel_->getLinearVelocity( );
    }

    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        motionModel_->setLinearVelocity( parameterValue );
    }

    int getParameterSize( )
    {
        return 3;
    }

private:
    std::shared_ptr< ground_stations::LinearStationMotionModel > motionModel_;
};

}  // namespace estimatable_parameters

}  // namespace tudat

#endif  // TUDAT_GROUNDSTATIONLINEARVELOCITY_H

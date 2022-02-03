/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_CUSTOM_ACCELERATION_MODEL_H
#define TUDAT_CUSTOM_ACCELERATION_MODEL_H

#include "tudat/astro/basic_astro/accelerationModel.h"

namespace tudat
{  
namespace basic_astrodynamics
{
class CustomAccelerationModel: public basic_astrodynamics::AccelerationModel3d
{
public:
    CustomAccelerationModel(
            const std::function< Eigen::Vector3d( const double ) > accelerationFunction ):
        accelerationFunction_( accelerationFunction )
    {
    }

    virtual void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            currentAcceleration_ = accelerationFunction_( currentTime );
        }

    }

private:
    std::function< Eigen::Vector3d( const double ) > accelerationFunction_;
};


} // namespace basic_astrodynamics
} // namespace tudat

#endif // TUDAT_CUSTOM_ACCELERATION_MODEL_H

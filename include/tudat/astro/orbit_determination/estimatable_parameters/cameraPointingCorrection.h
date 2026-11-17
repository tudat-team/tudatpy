/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CAMERA_POINTING_CORRECTION_H
#define TUDAT_CAMERA_POINTING_CORRECTION_H

#include <Eigen/Core>

#include "tudat/astro/system_models/vehicleSystems.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

class CameraPointingCorrection : public EstimatableParameter< Eigen::VectorXd >
{
public:
    CameraPointingCorrection( const std::shared_ptr< system_models::VehicleSystems > systemModels,
                              const std::string& associatedBody,
                              const std::string& cameraName ):
        EstimatableParameter< Eigen::VectorXd >( camera_pointing_correction, associatedBody, cameraName ), systemModels_( systemModels )
    {
        if( systemModels_ == nullptr )
        {
            throw std::runtime_error( "Error when making camera pointing correction parameter for " + associatedBody + ", " + cameraName +
                                      ": vehicle systems are null." );
        }
        if( systemModels_->getCameraMap( ).count( cameraName ) == 0 )
        {
            throw std::runtime_error( "Error when making camera pointing correction parameter for " + associatedBody + ", " + cameraName +
                                      ": camera not found." );
        }
    }

    ~CameraPointingCorrection( ) {}

    Eigen::VectorXd getParameterValue( )
    {
        return systemModels_->getCamera( parameterName_.second.second )->getPointingCorrection( );
    }

    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        systemModels_->getCamera( parameterName_.second.second )->setPointingCorrection( parameterValue );
    }

    //! Function to get the parameter size: three per picture of the camera.
    int getParameterSize( )
    {
        return 3 * systemModels_->getCamera( parameterName_.second.second )->getNumberOfPointingCorrections( );
    }

private:
    std::shared_ptr< system_models::VehicleSystems > systemModels_;
};

}  // namespace estimatable_parameters

}  // namespace tudat

#endif  // TUDAT_CAMERA_POINTING_CORRECTION_H

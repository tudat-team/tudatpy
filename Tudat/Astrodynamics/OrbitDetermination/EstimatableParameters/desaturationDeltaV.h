/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DESATURATIONDELTAV_H
#define TUDAT_DESATURATIONDELTAV_H

#include "Tudat/Astrodynamics/Propulsion/thrustAccelerationModel.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

class DesaturationDeltaV: public EstimatableParameter< Eigen::VectorXd >
{

public:


    DesaturationDeltaV(
            const std::shared_ptr< propulsion::MomentumWheelDesaturationThrustAcceleration > accelerationModel,
            const std::string& associatedBody ):
        EstimatableParameter< Eigen::VectorXd >( desaturation_delta_v_values, associatedBody ),
        accelerationModel_( accelerationModel )
    {
        numberOfDeltaVBlocks_  = accelerationModel_->getDeltaVValues( ).size( );
    }

    ~DesaturationDeltaV( ) { }

    Eigen::VectorXd getParameterValue( )
    {
        std::vector< Eigen::Vector3d > deltaVList = accelerationModel_->getDeltaVValues( );
        Eigen::VectorXd deltaVVector = Eigen::VectorXd( 3 * numberOfDeltaVBlocks_ );

        for( int i = 0; i < numberOfDeltaVBlocks_; i++ )
        {
           deltaVVector.segment( i * 3, 3 ) = deltaVList.at( i );
        }
        return deltaVVector;
    }

    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        std::vector< Eigen::Vector3d > deltaVList;
        for( int i = 0; i < numberOfDeltaVBlocks_; i++ )
        {
            deltaVList[ i ] = parameterValue.segment( i * 3, 3 );
        }

        accelerationModel_->setDeltaVValues( deltaVList );

    }

    int getParameterSize( )
    {
        return 3 * numberOfDeltaVBlocks_;
    }

protected:

private:

    const std::shared_ptr< propulsion::MomentumWheelDesaturationThrustAcceleration > accelerationModel_;

    int numberOfDeltaVBlocks_;

};

} // namespace estimatable_parameters

} // namespace tudat


#endif // TUDAT_DESATURATIONDELTAV_H

/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_AREATOMASSSCALINGFACTOR_H
#define TUDAT_AREATOMASSSCALINGFACTOR_H

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

class AreaToMassScalingFactor : public EstimatableParameter< double >
{
public:
    AreaToMassScalingFactor( const std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > > accelerationModels,
                             const std::string& associatedBody ):
        EstimatableParameter< double >( area_to_mass_scaling_factor, associatedBody ), accelerationModels_( accelerationModels )
    {
        double currentScalingFactor = accelerationModels_.at( 0 )->getAccelerationScalingFactor( );
        for( unsigned int i = 1; i < accelerationModels_.size( ); i++ )
        {
            if( accelerationModels_.at( i )->getAccelerationScalingFactor( ) != currentScalingFactor )
            {
                std::cerr << "Warning when creating area to mass scaling factor, scaling factors of constituent acceleration models are not initialized to the same value, setting all values to " +
                                std::to_string( currentScalingFactor )
                          << std::endl;
                accelerationModels_.at( i )->setAccelerationScalingFactor( currentScalingFactor );
            }
        }
    }

    ~AreaToMassScalingFactor( ) { }

    double getParameterValue( )
    {
        double currentScalingFactor = accelerationModels_.at( 0 )->getAccelerationScalingFactor( );
        for( unsigned int i = 1; i < accelerationModels_.size( ); i++ )
        {
            if( accelerationModels_.at( i )->getAccelerationScalingFactor( ) != currentScalingFactor )
            {
                throw std::runtime_error(
                        "Error when retrieving area to mass scaling factor value from parameters, scaling factors of constituent acceleration models is not consistent" );
            }
        }
        return currentScalingFactor;
    }

    void setParameterValue( double parameterValue )
    {
        for( unsigned int i = 0; i < accelerationModels_.size( ); i++ )
        {
            accelerationModels_.at( i )->setAccelerationScalingFactor( parameterValue );
        }
    }

    int getParameterSize( )
    {
        return 1;
    }

    std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > > getAccelerationModels( )
    {
        return accelerationModels_;
    }
protected:
    std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > > accelerationModels_;
};

class FullAccelerationScalingFactorParameter : public EstimatableParameter< double >
{
public:
    FullAccelerationScalingFactorParameter(
            const std::shared_ptr< basic_astrodynamics::AccelerationModel3d > accelerationModel,
            const std::string& bodyUndergoingAcceleration,
            const std::string& bodyExertingAcceleration ):
        EstimatableParameter< double >( full_acceleration_scaling_factor, bodyUndergoingAcceleration, bodyExertingAcceleration ),
        accelerationModel_( accelerationModel )
    {

    }

    ~FullAccelerationScalingFactorParameter( ) { }

    double getParameterValue( )
    {
        return accelerationModel_->getAccelerationScalingFactor(  );

    }

    void setParameterValue( double parameterValue )
    {
        accelerationModel_->setAccelerationScalingFactor( parameterValue );
    }

    int getParameterSize( )
    {
        return 1;
    }

    std::shared_ptr< basic_astrodynamics::AccelerationModel3d > getAccelerationModel( )
    {
        return accelerationModel_;
    }
protected:
    const std::shared_ptr< basic_astrodynamics::AccelerationModel3d > accelerationModel_;

};


}  // namespace estimatable_parameters

}  // namespace tudat

#endif  // TUDAT_AREATOMASSSCALINGFACTOR_H
/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_AERODYNAMICSCALINGCOEFFICIENT_H
#define TUDAT_AERODYNAMICSCALINGCOEFFICIENT_H

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

class AerodynamicScalingFactor : public EstimatableParameter< double >
{
public:
    AerodynamicScalingFactor( const std::shared_ptr< aerodynamics::AerodynamicAcceleration > aerodynamicAcceleration,
                              const EstimatebleParametersEnum parameterType,
                              const std::string& associatedBody ):
        EstimatableParameter< double >( parameterType, associatedBody ),
        aerodynamicAcceleration_( aerodynamicAcceleration )
    {
        if( ( parameterType != drag_component_scaling_factor ) &&
            ( parameterType != side_component_scaling_factor ) && 
            ( parameterType != lift_component_scaling_factor ) )
        {
            throw std::runtime_error( "Error when creating aerodynamic scaling parameter, type is inconsistent: " +
                                      std::to_string( parameterType ) );
        }
    }

    ~AerodynamicScalingFactor( ) { }

    double getParameterValue( )
    {
        if( parameterName_.first == drag_component_scaling_factor )
        {
            return aerodynamicAcceleration_->getDragComponentScaling( );
        }
        else if( parameterName_.first == side_component_scaling_factor )
        {
            return aerodynamicAcceleration_->getSideComponentScaling( );
        }
        else if( parameterName_.first == lift_component_scaling_factor )
        {
            return aerodynamicAcceleration_->getLiftComponentScaling( );
        }
        else
        {
            throw std::runtime_error( "Error when getting aerodynamic scaling parameter, type is inconsistent: " +
                                      std::to_string( parameterName_.first ) );
        }
    }

    void setParameterValue( double parameterValue )
    {
        if( parameterName_.first == drag_component_scaling_factor )
        {
            return aerodynamicAcceleration_->setDragComponentScaling( parameterValue );
        }
        else if( parameterName_.first == side_component_scaling_factor )
        {
            return aerodynamicAcceleration_->setSideComponentScaling( parameterValue );
        }
        else if( parameterName_.first == lift_component_scaling_factor )
        {
            return aerodynamicAcceleration_->setLiftComponentScaling( parameterValue );
        }
        else
        {
            throw std::runtime_error( "Error when setting aerodynamic scaling parameter, type is inconsistent: " +
                                      std::to_string( parameterName_.first ) );
        }
    }

    int getParameterSize( )
    {
        return 1;
    }

protected:
    std::shared_ptr< aerodynamics::AerodynamicAcceleration > aerodynamicAcceleration_;
};

}  // namespace estimatable_parameters

}  // namespace tudat

#endif  // TUDAT_AERODYNAMICSCALINGCOEFFICIENT_H
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
    AerodynamicScalingFactor( const std::vector< std::shared_ptr< aerodynamics::AerodynamicAcceleration > > aerodynamicAccelerations,
                              const EstimatebleParametersEnum parameterType,
                              const std::string& associatedBody ):
        EstimatableParameter< double >( parameterType, associatedBody ), aerodynamicAccelerations_( aerodynamicAccelerations )
    {
        if( ( parameterType != drag_component_scaling_factor ) && ( parameterType != side_component_scaling_factor ) &&
            ( parameterType != lift_component_scaling_factor ) )
        {
            throw std::runtime_error( "Error when creating aerodynamic scaling parameter, type is inconsistent: " +
                                      std::to_string( parameterType ) );
        }

        switch( parameterType )
        {
            case drag_component_scaling_factor:
                parameterIndex_ = 0;
                break;
            case side_component_scaling_factor:
                parameterIndex_ = 1;
                break;
            case lift_component_scaling_factor:
                parameterIndex_ = 2;
                break;
            default:
                throw std::runtime_error( "Error when creating aerodynamic scaling parameter, type is inconsistent: " +
                                          std::to_string( parameterType ) );
        }
    }

    AerodynamicScalingFactor( const std::shared_ptr< aerodynamics::AerodynamicAcceleration > aerodynamicAcceleration,
                              const EstimatebleParametersEnum parameterType,
                              const std::string& associatedBody ):
        AerodynamicScalingFactor( std::vector< std::shared_ptr< aerodynamics::AerodynamicAcceleration > >( { aerodynamicAcceleration } ),
                                  parameterType,
                                  associatedBody )
    {}

    ~AerodynamicScalingFactor( ) {}

    double getParameterValue( )
    {
        double parameterValue = aerodynamicAccelerations_.at( 0 )->getComponentScaling( parameterIndex_ );
        for( unsigned int i = 1; i < aerodynamicAccelerations_.size( ); i++ )
        {
            if( aerodynamicAccelerations_.at( i )->getComponentScaling( parameterIndex_ ) != parameterValue )
            {
                std::cerr << "Warning when retrieving aerodynamic acceleration scaling factor from list. List entries are not equal, "
                             "returning first entry"
                          << std::endl;
                break;
            }
        }
        return parameterValue;
    }

    void setParameterValue( double parameterValue )
    {
        for( unsigned int i = 0; i < aerodynamicAccelerations_.size( ); i++ )
        {
            aerodynamicAccelerations_.at( i )->setComponentScaling( parameterValue, parameterIndex_ );
        }
    }

    int getParameterSize( )
    {
        return 1;
    }

protected:
    std::vector< std::shared_ptr< aerodynamics::AerodynamicAcceleration > > aerodynamicAccelerations_;

    int parameterIndex_;
};

}  // namespace estimatable_parameters

}  // namespace tudat

#endif  // TUDAT_AERODYNAMICSCALINGCOEFFICIENT_H
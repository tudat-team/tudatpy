/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RELATIVEGRAVITATIONALPARAMETERRATE_H
#define TUDAT_RELATIVEGRAVITATIONALPARAMETERRATE_H

#include "tudat/astro/gravitation/gravitationalParameterRateAcceleration.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Relative rate (d mu / dt) / mu of the gravitational parameter of a body, in all gravitational parameter rate accelerations
//! that it exerts.
class RelativeGravitationalParameterRate : public EstimatableParameter< double >
{
public:
    RelativeGravitationalParameterRate(
            const std::vector< std::shared_ptr< gravitation::GravitationalParameterRateAcceleration > > accelerationModels,
            const std::string& associatedBody ):
        EstimatableParameter< double >( relative_gravitational_parameter_rate, associatedBody ), accelerationModels_( accelerationModels )
    {
        for( unsigned int i = 1; i < accelerationModels_.size( ); i++ )
        {
            if( accelerationModels_.at( i )->getRelativeGravitationalParameterRate( ) !=
                accelerationModels_.at( 0 )->getRelativeGravitationalParameterRate( ) )
            {
                throw std::runtime_error( "Error when making relative gravitational parameter rate parameter of " + associatedBody +
                                          ", rates are different in the acceleration models." );
            }
        }
    }

    double getParameterValue( )
    {
        return accelerationModels_.at( 0 )->getRelativeGravitationalParameterRate( );
    }

    void setParameterValue( double parameterValue )
    {
        for( unsigned int i = 0; i < accelerationModels_.size( ); i++ )
        {
            accelerationModels_.at( i )->setRelativeGravitationalParameterRate( parameterValue );
        }
    }

    int getParameterSize( )
    {
        return 1;
    }

private:
    std::vector< std::shared_ptr< gravitation::GravitationalParameterRateAcceleration > > accelerationModels_;
};

}  // namespace estimatable_parameters

}  // namespace tudat

#endif  // TUDAT_RELATIVEGRAVITATIONALPARAMETERRATE_H

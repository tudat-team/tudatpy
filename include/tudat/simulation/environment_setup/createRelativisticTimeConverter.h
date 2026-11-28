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

#ifndef TUDAT_CREATE_RELATIVISTIC_TIME_CONVERTER_H
#define TUDAT_CREATE_RELATIVISTIC_TIME_CONVERTER_H

#include <cmath>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/relativisticTimeConverterSettings.h"
#include "tudat/simulation/environment_setup/setRelativisticTimeConverters.h"
#include "tudat/simulation/propagation_setup/singleArcDynamicsSimulator.h"

namespace tudat
{

namespace simulation_setup
{

template< typename StateScalarType, typename TimeType >
void setRelativisticTimeConverter(
        const std::shared_ptr< DirectRelativisticTimeConverterSettings< StateScalarType, TimeType > >& conversionSettings,
        const SystemOfBodies& bodyMap )
{
    const auto barycentricSettings = conversionSettings->getBaryCentricToBodyCentricConversionSettings( );
    const auto topocentricConversions = conversionSettings->getBodyCentricToTopocentricConversionSettings( );

    const TimeType initialTime = barycentricSettings->getInitialTime( );
    if( !std::isfinite( static_cast< double >( initialTime ) ) )
    {
        throw std::runtime_error( "Error in setRelativisticTimeConverter: initial propagation time is not finite." );
    }

    const auto baseIntegratorSettings = conversionSettings->getNumericalIntegrationSettings( );
    if( baseIntegratorSettings == nullptr )
    {
        throw std::runtime_error( "Error in setRelativisticTimeConverter: no integrator settings provided." );
    }

    barycentricSettings->getOutputSettings( )->setIntegratedResult( true );

    auto barycentricIntegratorSettings = baseIntegratorSettings->clone( );
    barycentricSettings->setIntegratorSettings( barycentricIntegratorSettings );

    propagators::SingleArcDynamicsSimulator< StateScalarType, TimeType > barycentricSimulator( bodyMap, barycentricSettings, true );

    for( const auto& topocentricSettings : topocentricConversions )
    {
        if( topocentricSettings->getRelativisticStateDerivativeType( ) != propagators::first_order_bodycentric_to_topocentric )
        {
            throw std::runtime_error( "Error in setRelativisticTimeConverter: inconsistent derivative type for topocentric conversion." );
        }

        const TimeType topocentricInitialTime = topocentricSettings->getInitialTime( );
        if( !std::isfinite( static_cast< double >( topocentricInitialTime ) ) )
        {
            throw std::runtime_error( "Error in setRelativisticTimeConverter: topocentric propagation time is not finite." );
        }

        topocentricSettings->getOutputSettings( )->setIntegratedResult( true );

        auto topocentricIntegratorSettings = baseIntegratorSettings->clone( );
        topocentricSettings->setIntegratorSettings( topocentricIntegratorSettings );

        propagators::SingleArcDynamicsSimulator< StateScalarType, TimeType > topocentricSimulator( bodyMap, topocentricSettings, true );
    }
}

template< typename StateScalarType, typename TimeType >
void setRelativisticTimeConverters(
        const SystemOfBodies& bodyMap,
        const std::map< std::string, std::shared_ptr< DirectRelativisticTimeConverterSettings< StateScalarType, TimeType > > >&
                converterSettings )
{
    for( const auto& converterEntry : converterSettings )
    {
        setRelativisticTimeConverter( converterEntry.second, bodyMap );
    }
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_CREATE_RELATIVISTIC_TIME_CONVERTER_H

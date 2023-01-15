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
//RP-OLD

#include <iostream>

#include "tudat/astro/basic_astro/missionGeometry.h"
#include "tudat/astro/electromagnetism/radiationPressureInterface.h"


namespace tudat
{

namespace electromagnetism
{

//! Calculate radiation pressure at certain distance from a source.
double calculateRadiationPressure( const double sourcePower, const double distanceFromSource )
{
    return sourcePower / ( 4.0 * mathematical_constants::PI * distanceFromSource *
                           distanceFromSource * physical_constants::SPEED_OF_LIGHT );
}

void RadiationPressureInterface::updateInterfaceBase(
        const double currentTime )
{
    currentTime_ = currentTime;

    // Calculate current radiation pressure
    currentSolarVector_ = sourcePositionFunction_( ) - targetPositionFunction_( );
    double distanceFromSource = currentSolarVector_.norm( );
    currentRadiationPressure_ = calculateRadiationPressure(
                sourcePower_( ), distanceFromSource );

    // Calculate total shadowing due to occulting body; note that multiple concurrent
    // occultations are not completely correctly (prints warning).
    double shadowFunction = 1.0;
    double currentShadowFunction = 1.0;
    for( unsigned int i = 0; i < occultingBodyPositions_.size( ); i++ )
    {
        currentShadowFunction *= mission_geometry::computeShadowFunction(
                    sourcePositionFunction_( ), sourceRadius_, occultingBodyPositions_[ i ]( ),
                    occultingBodyRadii_[ i ], targetPositionFunction_( ) );

        if( currentShadowFunction != 1.0 && shadowFunction != 1.0 )
        {
            std::cerr << "Warning, multiple occultation occured in radiation pressure interface, results may be slightly in error" << std::endl;
        }

        shadowFunction *= currentShadowFunction;
    }

    currentRadiationPressure_ *= shadowFunction;
}

//! Function to update the current value of the radiation pressure
void RadiationPressureInterface::updateInterface(
        const double currentTime )
{   
    updateInterfaceBase( currentTime );
    radiationPressureCoefficient_ = radiationPressureCoefficientFunction_( currentTime );
}

} // namespace electromagnetism
} // namespace tudat

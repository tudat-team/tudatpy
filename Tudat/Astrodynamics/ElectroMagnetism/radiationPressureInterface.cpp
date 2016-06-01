/*    Copyright (c) 2010-2014, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      150408    D. Dirkx          File created.
 *
 *    References
 *
 *    Notes
 *
 */

#include <iostream>

#include "Tudat/Astrodynamics/BasicAstrodynamics/missionGeometry.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/radiationPressureInterface.h"


namespace tudat
{

namespace electro_magnetism
{

//! Calculate radiation pressure at certain distance from a source.
double calculateRadiationPressure( const double sourcePower, const double distanceFromSource )
{
    return sourcePower / ( 4.0 * mathematical_constants::PI * distanceFromSource *
                           distanceFromSource * physical_constants::SPEED_OF_LIGHT );
}

//! Function to update the current value of the radiation pressure
void RadiationPressureInterface::updateInterface(
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

} // namespace electro_magnetism
} // namespace tudat

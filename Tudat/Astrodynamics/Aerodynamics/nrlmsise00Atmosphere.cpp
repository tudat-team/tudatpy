/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *
 *    References
 *
 *    Notes
 *
 */
#include <iostream>
#include "Tudat/Astrodynamics/Aerodynamics/nrlmsise00Atmosphere.h"


//! Tudat library namespace.
namespace tudat
{
namespace aerodynamics
{

//! Compute the local atmospheric properties.
void NRLMSISE00Atmosphere::computeProperties(double altitude, double longitude,
                                             double latitude, double time) {
    // Compute the hash key
    size_t hashKey = hashFunc(altitude, longitude, latitude, time);

    // If hash key is same do nothing
    if (hashKey == hashKey_)
    {
      return;
    }
    hashKey_ = hashKey;

    // Retrieve input data.
    NRLMSISE00Input inputData = nrlmsise00InputFunction_(
        altitude, longitude, latitude, time );
    std::copy( inputData.apVector.begin( ), inputData.apVector.end( ), aph_.a );
    std::copy( inputData.switches.begin( ), inputData.switches.end(), flags_.switches);

    // Set environment properties.
    input_.g_lat = latitude;
    input_.g_long = longitude;
    input_.alt = altitude;
    input_.year = inputData.year;
    input_.doy = inputData.dayOfTheYear;
    input_.sec = inputData.secondOfTheDay;
    input_.lst = inputData.localSolarTime;
    input_.f107 = inputData.f107;
    input_.f107A = inputData.f107a;
    input_.ap = inputData.apDaily;
    input_.ap_a = &aph_;

    // Call Neutral Atmosphere Empircial Model from the surface to lower exosphere from NRLMSISE00.
    gtd7( &input_, &flags_, &output_ );

    // Compute current density and temperature.
    density_ = output_.d[ 5 ] * 1000.0;
    temperature_ = output_.t[ 1 ];
    if( useIdealGasLaw_ )
    {
        pressure_ = density_ * physical_constants::SPECIFIC_GAS_CONSTANT_AIR * temperature_;
    }
    else
    {
        pressure_ = TUDAT_NAN;
    }

}

//! Get the full model output
std::pair< std::vector< double >, std::vector< double > >
NRLMSISE00Atmosphere::getFullOutput( const double altitude, const double longitude,
                                     const double latitude, const double time )
{
    // Compute the properties
    computeProperties(altitude, longitude, latitude, time);
    std::pair< std::vector< double >, std::vector< double >> output;

    // Copy array members of struct to vectors on the pair.
    output.first = std::vector< double >(
                output_.d, output_.d + sizeof output_.d / sizeof output_.d[ 0 ] );
    output.second = std::vector< double >(
                output_.t, output_.t + sizeof output_.t / sizeof output_.t[ 0 ] );
    return output;
}

}  // namespace aerodynamics
}  // namespace tudat

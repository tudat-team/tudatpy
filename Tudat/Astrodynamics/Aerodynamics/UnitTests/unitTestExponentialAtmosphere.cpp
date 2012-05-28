/*    Copyright (c) 2010-2012, Delft University of Technology
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
 *      110316    F.M. Engelen      File created.
 *      110324    J. Melman         Time taken out of the equation.
 *      110427    F.M. Engelen      Changed input arguments for functions.
 *      110629    F.M. Engelen      Simplified unit test, removed dependency on other classes.
 *      110705    F.M. Engelen      Update to relative numerical errors.
 *      111128    B. Tong Minh      Added location-independent function test.
 *      111211    K. Kumar          Minor corrections to location-independent function test.
 *
 *    References
 *
 */

#include <cmath>
#include <iostream>
#include <limits>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/physicalConstants.h>

#include "Tudat/Astrodynamics/Aerodynamics/exponentialAtmosphere.h"
#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolation.h"

//! Test of implementation of the exponential atmosphere.
int main( )
{

    // Using declarations.
    using std::cerr;
    using std::endl;
    using std::fabs;

    // Declare test variable.
    bool isExponentialAtmosphereErroneous = false;

    // Summary of tests tests.
    // Test 1: Test set- and get-functions of constants.
    // Test 2: Test exponential atmosphere at sea level.
    // Test 3: Test exponential atmosphere at 10 km altitude.
    // Test 4: Test if the position-independent functions work.

    // Test 1: Test set- and get-functions of constants.
    // Initialize constants that need to be set.
    double constantTemperature = 288.16;
    double densityAtZeroAltitude = 1.225;
    double scaleHeight = 7.050e3;

    // Create an exponential atmosphere object.
    tudat::ExponentialAtmosphere exponentialAtmosphere;

    // Initialize the exponential atmosphere.
    exponentialAtmosphere.setConstantTemperature( constantTemperature );
    exponentialAtmosphere.setDensityAtZeroAltitude( densityAtZeroAltitude );
    exponentialAtmosphere.setScaleHeight( scaleHeight );
    exponentialAtmosphere.setSpecificGasConstant(
                tudat::physical_constants::SPECIFIC_GAS_CONSTANT_AIR );

    // Check if set and get functions work well.
    if ( fabs( ( exponentialAtmosphere.getConstantTemperature( )
                 - constantTemperature ) / constantTemperature )
         > std::numeric_limits< double >::epsilon( )
         || fabs( ( exponentialAtmosphere.getDensityAtZeroAltitude( )
                    - densityAtZeroAltitude ) / densityAtZeroAltitude )
         > std::numeric_limits< double >::epsilon( )
         || fabs( ( exponentialAtmosphere.getScaleHeight( ) - scaleHeight ) / scaleHeight )
         > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The get or set functions for the constants of the exponential atmosphere ";
        cerr << "do not work correctly." << endl;
        isExponentialAtmosphereErroneous = true;
    }

    // Test 2: Check if the atmosphere is calculated correctly at sea level.
    // Values from "US Standard Atmosphere 1976,
    // http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770009539_1977009539.pdf".

    // Initialize the altitude.
    double altitude = 0.0;

    // Initialize the pressure at zero altitude.
    double pressureAtZeroAltitude = 101325.0;

    // Check whether the atmosphere is calculated correctly at sea level.
    if ( fabs( ( exponentialAtmosphere.getTemperature( altitude )
                 - constantTemperature ) / constantTemperature )
         > std::numeric_limits< double >::epsilon( )
         || fabs( ( exponentialAtmosphere.getDensity( altitude )
                    - densityAtZeroAltitude ) / densityAtZeroAltitude )
         > std::numeric_limits< double >::epsilon( )
         || fabs( ( exponentialAtmosphere.getPressure( altitude )
                    - pressureAtZeroAltitude ) / pressureAtZeroAltitude )
         > 0.002 * pressureAtZeroAltitude )
        // Because of different gas constant used in the USSA1976, there is a slight difference.
    {
        cerr << "The exponential atmosphere at sea level is calculated incorrectly." << endl;
        isExponentialAtmosphereErroneous = true;
    }

    // Test 3: Test exponential atmosphere at 10 km altitude.
    // Check whether the atmosphere is calculated correctly at 10 km.
    // The given value for pressure was calculated off-line with a calculator.
    altitude = 10.0e3;

    // Also use longitude and latitude to check overloading.
    double longitude = 0.0;
    double latitude = 0.0;
    double time = 0.0;

    // Declare and set expected density.
    double expectedDensity  = densityAtZeroAltitude * std::exp ( -altitude / scaleHeight );

    if ( fabs( ( exponentialAtmosphere.getTemperature( altitude, longitude, latitude, time )
                 - constantTemperature ) / constantTemperature )
         > std::numeric_limits< double >::epsilon( )
         || fabs( ( exponentialAtmosphere.getDensity( altitude, longitude, latitude, time )
                    - expectedDensity ) / expectedDensity )
         > std::numeric_limits< double >::epsilon( )
         || fabs( ( exponentialAtmosphere.getPressure( altitude, longitude, latitude, time )
                    - 24526.24934607106 ) ) > 1.0e-10 )
    {
        cerr << "The exponential atmosphere at 10 km altitude is calculated incorrectly." << endl;
        isExponentialAtmosphereErroneous = true;
    }

    // Test 4: Test if the position-independent functions work.
    double density1 = exponentialAtmosphere.getDensity( altitude );
    double density2 = exponentialAtmosphere.getDensity( altitude, longitude, latitude, time );

    double pressure1 = exponentialAtmosphere.getPressure( altitude );
    double pressure2 = exponentialAtmosphere.getPressure( altitude, longitude, latitude, time );

    double temperature1 = exponentialAtmosphere.getTemperature( altitude );
    double temperature2 = exponentialAtmosphere.getTemperature( altitude, longitude, latitude, time );

    if ( fabs( ( density1 - density2 ) / density2 ) > std::numeric_limits< double >::epsilon( )
         || fabs( ( pressure1 - pressure2 ) / pressure2 ) > std::numeric_limits< double >::epsilon( )
         || fabs( ( temperature1 - temperature2 ) /temperature2 ) > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "Location-dependent and location-independent functions did not give the same "
             << "result." << endl;

        // For some reason if you don't use the temporary variables, the optimizer will do weird
        // stuff causing density1 != density2 to hold true...

        cerr << "Density difference: " << ( density1 - density2 ) << endl
             << "Pressure difference: " << ( pressure1 - pressure2 ) << endl
             << "Temperature difference: " << ( temperature1 - temperature2 ) << endl;

        isExponentialAtmosphereErroneous = true;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    if ( isExponentialAtmosphereErroneous )
    {
        cerr << "testExponentialAtmosphere failed!" << std::endl;
    }

    return isExponentialAtmosphereErroneous;
}

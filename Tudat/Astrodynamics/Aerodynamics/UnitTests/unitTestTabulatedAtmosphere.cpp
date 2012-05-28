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
 *      110713    F.M. Engelen      File created.
 *      110721    J. Melman         Alignment, comments, error messages, and consistency modified.
 *      110722    F.M. Engelen      Replaced values in to book values.
 *      111128    B. Tong Minh      Added location-independent function test.
 *      111211    K. Kumar          Minor corrections to location-independent function test.
 *
 *    References
 *      Introduction to Flight, Fifth edition, Appendix A, John D. Anderson Jr., McGraw Hill, 2005.
 *
 */

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "Tudat/Astrodynamics/Aerodynamics/tabulatedAtmosphere.h"
#include "Tudat/InputOutput/basicInputOutput.h"

//! Test implementation of the tabulated atmosphere.
int main( )
{
    // Using declarations.
    using std::cerr;
    using std::endl;
    using std::pow;
    using std::fabs;
    using tudat::TabulatedAtmosphere;

    // Declare test variable.
    bool isTabulatedAtmosphereErroneous = false;

    // Test 1: Test tabulated atmosphere at sea level.
    // Test 2: Test tabulated atmosphere at 10.0 km altitude including passing arbitrary longitude
    //         and latitude.
    // Test 3: Test tabulated atmosphere at 10.05 km altitude when just passing the altitude..
    // Test 4: Test tabulated atmosphere at 1000 km altitude with figure.
    // Test 5: Test tabulated atmosphere at 1000 km altitude with table.
    // Test 6: Test if the atmosphere file can be read multiple times.
    // Test 7: Test if the position-independent functions work.

    // Create a tabulated atmosphere object.
    TabulatedAtmosphere tabulatedAtmosphere;

    // Initialize it with the desired file.
    tabulatedAtmosphere.initialize( tudat::input_output::getTudatRootPath( ) +
                                    "Astrodynamics/Aerodynamics/AtmosphereTables/" +
                                    "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat" );

    // Test 1: Check if the atmosphere is calculated correctly at sea level.
    // Values from "US Standard Atmosphere 1976,
    // http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770009539_1977009539.pdf".

    // Initialize the altitude.
    double altitude = 0.0;

    // Check whether the atmosphere is calculated correctly at sea level.
    if ( fabs( ( tabulatedAtmosphere.getTemperature( altitude ) - 288.15 ) / 288.15 )
         > std::numeric_limits< double >::epsilon( )
         || fabs( ( tabulatedAtmosphere.getDensity( altitude ) - 1.225 ) / 1.225 )
         > std::numeric_limits< double >::epsilon( )
         || fabs( ( tabulatedAtmosphere.getPressure( altitude ) - 101325.0 ) / 101325.0 )
         > 1.0e-4 )
        // Because of different gas constant used in the USSA1976, there is a slight difference.
    {
        cerr << "The tabulated atmosphere at sea level is calculated incorrectly." << endl;
        cerr << "Temperature = " << tabulatedAtmosphere.getTemperature( altitude ) << " K"<< endl;
        cerr << "Density = " << tabulatedAtmosphere.getDensity( altitude ) << " kg/m^3" << endl;
        cerr << "Pressure = " << tabulatedAtmosphere.getPressure( altitude ) << " N/m^2" << endl;

        isTabulatedAtmosphereErroneous = true;
    }

    // Test 2: Test tabulated atmosphere at 10.0 km altitude including passing arbitrary longitude
    //         and latitude.
    // Check whether the atmosphere is calculated correctly at 10 km.
    // The given value for pressure was obtained from table in book
    altitude = 10.0e3 ;

    // Also use longitude and latitude to check overloading.
    double longitude = 0.0;
    double latitude = 0.0;
    double time = 0.0;

    if ( fabs( tabulatedAtmosphere.getTemperature( altitude, longitude, latitude, time ) - 223.26 )
         > 1.0e-2
         || fabs( tabulatedAtmosphere.getDensity( altitude, longitude, latitude, time ) - 0.41351 )
         > 1.0e-4
         || fabs( tabulatedAtmosphere.getPressure( altitude, longitude, latitude, time ) - 26500  )
         > 1.0 )
    {
        cerr << "The tabulated atmosphere at 10 km altitude is calculated incorrectly." << endl;
        cerr << "Temperature = " << tabulatedAtmosphere.getTemperature( altitude ) << " K."<< endl;
        cerr << "Density = " << tabulatedAtmosphere.getDensity( altitude ) << " kg/m^3." << endl;
        cerr << "Pressure = " << tabulatedAtmosphere.getPressure( altitude ) << " N/m^2." << endl;

        isTabulatedAtmosphereErroneous = true;
    }

    // Test 3: Test tabulated atmosphere at 10.05 km altitude when just passing the altitude..
    // Check whether the atmosphere is calculated correctly at 10.05 km.
    // The values are linear interpolated values based on book values.
    altitude = 10.05e3;

    if ( fabs(  tabulatedAtmosphere.getTemperature( altitude ) - 222.9350 ) > 2.0e-2
         || fabs(  tabulatedAtmosphere.getDensity( altitude ) - 0.4110 ) > 1.0e-3
         || fabs(  tabulatedAtmosphere.getPressure( altitude ) - 26299 ) > 1.0 )
    {
        cerr << "The tabulated atmosphere at 10.05 km altitude is calculated incorrectly." << endl;
        cerr << "Temperature = " << tabulatedAtmosphere.getTemperature( altitude ) << " K"<< endl;
        cerr << "Density = " << tabulatedAtmosphere.getDensity( altitude ) << " kg/m^3" << endl;
        cerr << "Pressure = " << tabulatedAtmosphere.getPressure( altitude ) << " N/m^2" << endl;

        isTabulatedAtmosphereErroneous = true;
    }

    // Test 4: Test tabulated atmosphere at 1000 km altitude.
    // Check whether the atmosphere is calculated correctly at 1000 km.
    // Compared with a figure in the US 1976 standard atmosphere document. (coars figure)
    altitude = 1.0e6 ;

    if ( fabs( ( tabulatedAtmosphere.getTemperature( altitude ) - 1000.0 ) / 1000.0 )
         > std::numeric_limits< double >::epsilon( )
         || fabs( ( tabulatedAtmosphere.getDensity( altitude ) - 5.0e-15 ) / 5.0e-15 ) > 0.5
         || fabs( ( tabulatedAtmosphere.getPressure( altitude ) - 1.0e-8 ) / 1.0e-8  ) > 0.5 )
    {
        cerr << "The tabulated atmosphere at 1000 km altitude is calculated incorrectly." << endl;
        cerr << "Temperature = " << tabulatedAtmosphere.getTemperature( altitude ) << " K"<< endl;
        cerr << "Density = " << tabulatedAtmosphere.getDensity( altitude ) << " kg/m^3" << endl;
        cerr << "Pressure = " << tabulatedAtmosphere.getPressure( altitude ) << " N/m^2" << endl;
        isTabulatedAtmosphereErroneous = true;
    }

    // Test 5: Test tabulated atmosphere at 1000 km altitude.
    // Check whether the atmosphere is calculated correctly at 1000 km.
    // Compared with the input table.
    altitude = 1.0e6;

    if ( fabs( ( tabulatedAtmosphere.getTemperature( altitude ) - 1000.0 ) / 1000.0 )
         > std::numeric_limits< double >::epsilon( )
         || fabs( ( tabulatedAtmosphere.getDensity( altitude ) - 3.5618e-15 ) / 3.5618e-15 )
         > std::numeric_limits< double >::epsilon( )
         || fabs( ( tabulatedAtmosphere.getPressure( altitude ) - 7.5158e-9 ) / 7.5158e-9 )
         > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The tabulated atmosphere at 1000 km altitude is calculated incorrectly." << endl;
        cerr << "Temperature = " << tabulatedAtmosphere.getTemperature( altitude ) << " K"<< endl;
        cerr << "Density = " << tabulatedAtmosphere.getDensity( altitude ) << " kg/m^3" << endl;
        cerr << "Pressure = " << tabulatedAtmosphere.getPressure( altitude ) << " N/m^2" << endl;
        isTabulatedAtmosphereErroneous = true;
    }

    // Test 6: Test if the atmosphere file can be read multiple times.
    try
    {
        tabulatedAtmosphere.initialize( tudat::input_output::getTudatRootPath( ) +
                                        "Astrodynamics/Aerodynamics/AtmosphereTables/" +
                                        "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat" );
    }

    catch ( std::runtime_error multipleFileReadError )
    {
        cerr << "Caught exception while opening the atmosphere file for a second time." << endl;
        isTabulatedAtmosphereErroneous = true;
    }

    // Test 7: Test if the position-independent functions work.
    double density1 = tabulatedAtmosphere.getDensity( altitude );
    double density2 = tabulatedAtmosphere.getDensity( altitude, longitude, latitude, time );

    double pressure1 = tabulatedAtmosphere.getPressure( altitude );
    double pressure2 = tabulatedAtmosphere.getPressure( altitude, longitude, latitude, time );

    double temperature1 = tabulatedAtmosphere.getTemperature( altitude );
    double temperature2 = tabulatedAtmosphere.getTemperature( altitude, longitude,
                                                              latitude, time );

    if ( fabs( density1 - density2 ) > std::numeric_limits< double >::epsilon( )
         || fabs( pressure1 - pressure2 ) > std::numeric_limits< double >::epsilon( )
         || fabs( temperature1 - temperature2 ) > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "Location-dependent and location-independent functions did not give the same "
             << "result." << endl;

        // For some reason if you don't use the temporary variables, the optimizer will do weird
        // stuff causing density1 != density2 to hold true...

        cerr << "Density difference: " << ( density1 - density2 ) << endl
             << "Pressure difference: " << ( pressure1 - pressure2 ) << endl
             << "Temperature difference: " << ( temperature1 - temperature2 ) << endl;

        isTabulatedAtmosphereErroneous = true;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    if ( isTabulatedAtmosphereErroneous )
    {
        cerr << "testTabulatedAtmosphere failed!" << std::endl;
    }

    return isTabulatedAtmosphereErroneous;
}

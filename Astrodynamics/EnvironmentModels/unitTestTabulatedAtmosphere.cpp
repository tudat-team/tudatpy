/*! \file unitTestTabulatedAtmosphere.cpp
 *    Source file that defines the Tabulated atmosphere unit test included in Tudat.
 *
 *    Path              : /Astrodynamics/EnvironmentModels/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 13 July, 2011
 *    Last modified     : 22 July, 2011
 *
 *    References
 *    Introduction to Flight, Fifth edition, Appendix A, John D. Anderson Jr., McGraw Hill,
 *    2005
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110713    F.M. Engelen      File created.
 *      110721    J. Melman         Alignment, comments, error messages,
 *                                  and consistency modified.
 *      110722    F.M. Engelen      Replaced values in to book values.
 */

// Include statements.
#include <cmath>
#include <iostream>
#include "Astrodynamics/EnvironmentModels/tabulatedAtmosphere.h"
#include "Astrodynamics/EnvironmentModels/unitTestTabulatedAtmosphere.h"
#include "Mathematics/basicMathematicsFunctions.h"

//! Namespace for all unit tests.
namespace unit_tests
{

//! Test of implementation of the Tabulated atmosphere.
bool testTabulatedAtmosphere( )
{
    // Using declarations.
    using std::cerr;
    using std::endl;
    using std::pow;
    using std::fabs;
    using mathematics::MACHINE_PRECISION_DOUBLES;

    // Declare test variable.
    bool isTabulatedAtmosphereBad_ = false;

    // Test 1: Test tabulated atmosphere at sea level.
    // Test 2: Test tabulated atmosphere at 10.0 km altitude including
    //         passing arbitrary longitude and latitude.
    // Test 3: Test tabulated atmosphere at 10.1 km altitude when just
    //         passing the altitude.
    // Test 4: Test tabulated atmosphere at 1000 km altitude with figure.
    // Test 5: Test tabulated atmosphere at 1000 km altitude with table.

    // Create a tabulated atmosphere object.
    TabulatedAtmosphere tabulatedAtmosphere;

    // Initialize it with the desired file.
    tabulatedAtmosphere.initialize( "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat" );

    // Test 1: Check if the atmosphere is calculated correctly at sea level.
    // Values from "US Standard Atmosphere 1976,
    // http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770009539_1977009539.pdf".

    // Initialize the altitude.
    double altitude = 0.0;

    // Check whether the atmosphere is calculated correctly at sea level.
    if ( fabs( ( tabulatedAtmosphere.getTemperature( altitude ) - 288.15 ) / 288.15 )
         > MACHINE_PRECISION_DOUBLES
         || fabs( ( tabulatedAtmosphere.getDensity( altitude ) - 1.225 ) / 1.225 )
         > MACHINE_PRECISION_DOUBLES
         || fabs( ( tabulatedAtmosphere.getPressure( altitude ) - 101325.0 ) / 101325.0 )
         > 1.0e-4 )
        // Because of different gas constant used in the USSA1976, there is a slight difference.
    {
        cerr << "The tabulated atmosphere at sea level is calculated incorrectly." << endl;
        cerr << "Temperature = " << tabulatedAtmosphere.getTemperature( altitude ) << " K"<< endl;
        cerr << "Density = " << tabulatedAtmosphere.getDensity( altitude ) << " kg/m^3" << endl;
        cerr << "Pressure = " << tabulatedAtmosphere.getPressure( altitude ) << " N/m^2" << endl;

        isTabulatedAtmosphereBad_ = true;
    }

    // Test 2: Test tabulated atmosphere at 10 km altitude including
    //         passing arbitrary longitude and latitude.
    // Check whether the atmosphere is calculated correctly at 10 km.
    // The given value for pressure was obtained from table in book
    altitude = 10.0e3 ;

    // Also use longitude and latitude to check overloading.
    double longitude = 0.0;
    double latitude = 0.0;

    if ( fabs(  tabulatedAtmosphere.getTemperature( altitude, longitude, latitude ) - 223.26 )
         > 1.0e-2
         || fabs(  tabulatedAtmosphere.getDensity( altitude, longitude, latitude ) - 0.41351 )
         > 1.0e-4
         || fabs(  tabulatedAtmosphere.getPressure( altitude, longitude, latitude ) - 26500  )
         > 1.0 )
    {
        cerr << "The tabulated atmosphere at 10 km altitude is calculated incorrectly." << endl;
        cerr << "Temperature = " << tabulatedAtmosphere.getTemperature( altitude ) << " K"<< endl;
        cerr << "Density = " << tabulatedAtmosphere.getDensity( altitude ) << " kg/m^3" << endl;
        cerr << "Pressure = " << tabulatedAtmosphere.getPressure( altitude ) << " N/m^2" << endl;

        isTabulatedAtmosphereBad_ = true;
    }

    // Test 3: Test tabulated atmosphere at 10.05 km altitude when just
    //         passing the altitude..
    // Check whether the atmosphere is calculated correctly at 10.05 km.
    // The values are linear interpolated values based on book values.
    altitude = 10.05e3;

    if ( fabs(  tabulatedAtmosphere.getTemperature( altitude ) - 222.9350 )  > 2.0e-2
         || fabs(  tabulatedAtmosphere.getDensity( altitude ) - 0.4110 ) > 1.0e-3
         || fabs(  tabulatedAtmosphere.getPressure( altitude ) - 26299 ) > 1.0 )
    {
        cerr << "The tabulated atmosphere at 10.05 km altitude is calculated incorrectly." << endl;
        cerr << "Temperature = " << tabulatedAtmosphere.getTemperature( altitude ) << " K"<< endl;
        cerr << "Density = " << tabulatedAtmosphere.getDensity( altitude ) << " kg/m^3" << endl;
        cerr << "Pressure = " << tabulatedAtmosphere.getPressure( altitude ) << " N/m^2" << endl;

        isTabulatedAtmosphereBad_ = true;
    }

    // Test 4: Test tabulated atmosphere at 1000 km altitude.
    // Check whether the atmosphere is calculated correctly at 1000 km.
    // Compared with a figure in the US 1976 standard atmosphere document. (coars figure)
    altitude = 1.0e6 ;

    if ( fabs( ( tabulatedAtmosphere.getTemperature( altitude ) - 1000.0 ) / 1000.0 )
         > MACHINE_PRECISION_DOUBLES
         || fabs( ( tabulatedAtmosphere.getDensity( altitude ) - 5.0e-15 ) / 5.0e-15 ) > 0.5
         || fabs( ( tabulatedAtmosphere.getPressure( altitude ) - 1.0e-8 ) / 1.0e-8  ) > 0.5 )
    {
        cerr << "The tabulated atmosphere at 1000 km altitude is calculated incorrectly." << endl;
        cerr << "Temperature = " << tabulatedAtmosphere.getTemperature( altitude ) << " K"<< endl;
        cerr << "Density = " << tabulatedAtmosphere.getDensity( altitude ) << " kg/m^3" << endl;
        cerr << "Pressure = " << tabulatedAtmosphere.getPressure( altitude ) << " N/m^2" << endl;
        isTabulatedAtmosphereBad_ = true;
    }

    // Test 5: Test tabulated atmosphere at 1000 km altitude.
    // Check whether the atmosphere is calculated correctly at 1000 km.
    // Compared with the input table.
    altitude = 1.0e6;

    if ( fabs( ( tabulatedAtmosphere.getTemperature( altitude ) - 1000.0 ) / 1000.0 )
         > MACHINE_PRECISION_DOUBLES
         || fabs( ( tabulatedAtmosphere.getDensity( altitude ) - 3.5618e-15 ) / 3.5618e-15 )
         > MACHINE_PRECISION_DOUBLES
         || fabs( ( tabulatedAtmosphere.getPressure( altitude ) - 7.5158e-9 ) / 7.5158e-9 )
         > MACHINE_PRECISION_DOUBLES )
    {
        cerr << "The tabulated atmosphere at 1000 km altitude is calculated incorrectly." << endl;
        cerr << "Temperature = " << tabulatedAtmosphere.getTemperature( altitude ) << " K"<< endl;
        cerr << "Density = " << tabulatedAtmosphere.getDensity( altitude ) << " kg/m^3" << endl;
        cerr << "Pressure = " << tabulatedAtmosphere.getPressure( altitude ) << " N/m^2" << endl;
        isTabulatedAtmosphereBad_ = true;
    }

    return isTabulatedAtmosphereBad_;
}

}

// End of file.

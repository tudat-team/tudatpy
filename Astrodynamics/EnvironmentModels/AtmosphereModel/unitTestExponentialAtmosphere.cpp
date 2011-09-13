/*! \file unitTestExponentialAtmosphere.cpp
 *    Source file that defines the exponential atmosphere unit test included in Tudat.
 *
 *    Path              : /Astrodynamics/EnvironmentModels/AtmosphereModel/
 *    Version           : 2
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
 *    Date created      : 16 March, 2011
 *    Last modified     : 29 June, 2011
 *
 *    References
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
 *      110316    F.M. Engelen      File created.
 *      110324    J. Melman         Time taken out of the equation.
 *      110427    F.M. Engelen      Changed input arguments for functions.
 *      110629    F.M. Engelen      Simplified unit test, removed dependency on other classes.
 *      110705    F.M. Engelen      Update to relative numerical errors.
 */

// Include statements.
#include <cmath>
#include "Astrodynamics/physicalConstants.h"
#include "Astrodynamics/EnvironmentModels/AtmosphereModel/exponentialAtmosphere.h"
#include "Astrodynamics/EnvironmentModels/AtmosphereModel/unitTestExponentialAtmosphere.h"
#include "Mathematics/basicMathematicsFunctions.h"

//! Namespace for all unit tests.
namespace unit_tests
{

//! Test of implementation of the exponential atmosphere.
bool testExponentialAtmosphere( )
{
    // Using declarations.
    using std::cerr;
    using std::endl;
    using mathematics::computeAbsoluteValue;
    using mathematics::MACHINE_PRECISION_DOUBLES;

    // Declare test variable.
    bool isExponentialAtmosphereBad_ = false;

    // Three tests.
    // Test 1: Test set and get functions of constants.
    // Test 2: Test exponential atmosphere at sea level.
    // Test 3: Test exponential atmosphere at 10 km altitude.

    // Test 1: check set and get functions of constants.

    // Initialize constants that need to be set.
    double constantTemperature = 288.16;
    double densityAtZeroAltitude = 1.225;
    double scaleHeight = 7.050e3;

    // Create an exponential atmosphere object.
    ExponentialAtmosphere exponentialAtmosphere;

    // Initialize the exponential atmosphere.
    exponentialAtmosphere.setConstantTemperature( constantTemperature );
    exponentialAtmosphere.setDensityAtZeroAltitude( densityAtZeroAltitude );
    exponentialAtmosphere.setScaleHeight( scaleHeight );
    exponentialAtmosphere.setSpecificGasConstant( PhysicalConstants::SPECIFIC_GAS_CONSTANT_AIR );

    // Check if set and get functions work well.
    if ( computeAbsoluteValue( ( exponentialAtmosphere.getConstantTemperature( )
                               - constantTemperature ) / constantTemperature )
        > MACHINE_PRECISION_DOUBLES
        ||
        computeAbsoluteValue( ( exponentialAtmosphere.getDensityAtZeroAltitude( )
                              - densityAtZeroAltitude ) / densityAtZeroAltitude )
        > MACHINE_PRECISION_DOUBLES
        ||
        computeAbsoluteValue( ( exponentialAtmosphere.getScaleHeight( )  - scaleHeight ) /
                              scaleHeight )
        > MACHINE_PRECISION_DOUBLES )
    {
        cerr << "The get or set functions for the constants of the exponential atmosphere ";
        cerr << "do not work correctly." << endl;
        isExponentialAtmosphereBad_ = true;
    }

    // Test 2: Check if the atmosphere is calculated correctly at sea level.
    // Values from "US Standard Atmosphere 1976,
    // http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770009539_1977009539.pdf".

    // Initialize the altitude.
    double altitude = 0.0;

    // Initialize the pressure at zero altitude.
    double pressureAtZeroAltitude = 101325.0;

    // Check whether the atmosphere is calculated correctly at sea level.
    if ( computeAbsoluteValue( ( exponentialAtmosphere.getTemperature( altitude )
                               - constantTemperature ) / constantTemperature )
        > MACHINE_PRECISION_DOUBLES
        ||
        computeAbsoluteValue( ( exponentialAtmosphere.getDensity( altitude )
                              - densityAtZeroAltitude ) / densityAtZeroAltitude )
        > MACHINE_PRECISION_DOUBLES
        ||
        computeAbsoluteValue( ( exponentialAtmosphere.getPressure( altitude )
                              - pressureAtZeroAltitude ) / pressureAtZeroAltitude )
        > 0.002 * pressureAtZeroAltitude )
        // Because of different gas constant used in the USSA1976, there is a slight difference.
    {
        cerr << "The exponential atmosphere at sea level is calculated incorrectly." << endl;
        isExponentialAtmosphereBad_ = true;
    }

    // Test 3: Test exponential atmosphere at 10 km altitude.
    // Check whether the atmosphere is calculated correctly at 10 km.
    // The given value for pressure was calculated off-line with a calculator.
    altitude = 10.0e3 ;

    // Also use longitude and latitude to check overloading.
    double longitude = 0.0;
    double latitude = 0.0;

    double expectedDensity  = densityAtZeroAltitude * exp ( - altitude / scaleHeight );

    if ( computeAbsoluteValue( ( exponentialAtmosphere.getTemperature(
            altitude, longitude, latitude )
                               - constantTemperature ) / constantTemperature )
        > MACHINE_PRECISION_DOUBLES
        ||
        computeAbsoluteValue( ( exponentialAtmosphere.getDensity(
                altitude, longitude, latitude )
                              - expectedDensity ) / expectedDensity )
        > MACHINE_PRECISION_DOUBLES
        ||
        computeAbsoluteValue( ( exponentialAtmosphere.getPressure(
                altitude, longitude, latitude )
                              - 24526.24934607106 ) )
        > 1.0e-10 )
    {
        cerr << "The exponential atmosphere at 10 km altitude is calculated incorrectly." << endl;
        isExponentialAtmosphereBad_ = true;
    }

    return isExponentialAtmosphereBad_;
}

}

// End of file.

/*! \file aerodynamics.cpp
 *  This file contains the definition of the aerodynamics namespace.
 *
 *  Path              : Astrodynamics/EnvironmentModels/Aerodynamics
 *  Version           : 3
 *  Check status      : Checked
 *
 *  Author            : Dominic Dirkx
 *  Affiliation       : TU Delft
 *  E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *  Checker           : L. Abdulkadir
 *  Affiliation       : TU Delft
 *  E-mail address    : L.Abdulkadir@student.tudelft.nl
 *
 *  Date created      : 25 November, 2010
 *  Last modified     : 10 February, 2011
 *
 *  References
 *
 *  Anderson Jr., J.D. , Fundamentals of Aerodynamics, 3rd edition, McGraw Hill
 *  2001.
 *  Gentry, A., Smyth, D., and Oliver, W. . The Mark IV Supersonic-Hypersonic
 *  Arbitrary Body Program, Volume II - Program Formulation, Douglas Aircraft
 *  Company, 1973.
 *  Anderson Jr. , J.D, Hypersonic and High-Temperature Gas Dynamics, 2nd
 *  edition, AIAA Education Series, 2006
 *
 *  Notes
 *
 *
 *  Copyright (c) 2010 Delft University of Technology.
 *
 *  This software is protected by national and international copyright.
 *  Any unauthorized use, reproduction or modification is unlawful and
 *  will be prosecuted. Commercial and non-private application of the
 *  software in any form is strictly prohibited unless otherwise granted
 *  by the authors.
 *
 *  The code is provided without any warranty; without even the implied
 *  warranty of merchantibility or fitness for a particular purpose.
 *
 *  Changelog
 *      YYMMDD    author        comment
 *      102511   D. Dirkx       First version of file.
 *      110501   D. Dirkx       Added more comments.
 *      110203   L. Abdulkadir  Code check.
 *      110208   D. Dirkx       Fixed shock temperature ratio.
 *      110210   L. Abdulkadir  Code check.
 */

#include "aerodynamics.h"

namespace aerodynamics{

//! Function to calculate local to static pressure ratio.
double calculateLocalToStaticPressureRatio(
        const double& machNumber,
        const double& ratioOfSpecificHeats)
{
    // Declare local variable.
    double localToStaticPressureRatio;

    // Calculate pressure ratio.
    localToStaticPressureRatio = pow( 2 / ( 2 + ( ratioOfSpecificHeats - 1 ) *
                            pow( machNumber, 2.0 ) ), ratioOfSpecificHeats /
                            ( ratioOfSpecificHeats - 1 ) );
    return localToStaticPressureRatio;
}

//! Function to evaluate the Prandtl-Meyer function.
double  evaluatePrandtlMeyerFunction( const double& machNumber,
                                      const double& ratioOfSpecificHeats )
{
    // Declare local variable.
    double prandtlMeyerFunction;

    // Evaluate Prandtl-Meyer function.
    prandtlMeyerFunction = sqrt ( ( ratioOfSpecificHeats + 1 ) /
        ( ratioOfSpecificHeats - 1 ) ) * atan ( sqrt ( ( ratioOfSpecificHeats - 1 ) /
        ( ratioOfSpecificHeats + 1 ) * ( pow( machNumber, 2.0 ) - 1 ) ) ) -
        atan( sqrt ( pow ( machNumber, 2.0 ) - 1 ) );
    return prandtlMeyerFunction;
}

//! Function to calculate stagnation pressure in supersonic flow, assuming flow
//! passes through normal shock.
double  calculateStagnationPressure(
        const double& machNumber,
        const double& ratioOfSpecificHeats)
{
    // Declare local variable.
    double stagnationPressureCoefficient ;

    // Calculate stagnation pressure coefficient.
    stagnationPressureCoefficient  = 2 / ( ratioOfSpecificHeats *
               pow( machNumber, 2.0 ) ) * ( pow ( pow ( ( ratioOfSpecificHeats + 1 ) *
               machNumber, 2.0 ) / ( 4 * ratioOfSpecificHeats * pow( machNumber , 2.0 ) -
               2 * ( ratioOfSpecificHeats - 1 ) ), ratioOfSpecificHeats /
              ( ratioOfSpecificHeats - 1 ) ) * ( ( 1 - ratioOfSpecificHeats +
              2 * ratioOfSpecificHeats * pow( machNumber , 2.0 ) ) /
                                                 ( ratioOfSpecificHeats + 1 ) ) - 1 );
    return stagnationPressureCoefficient;
}

//! Function to calculate pressure coefficient based on Newtonian theory.
double  calculateNewtonianPressureCoefficient(
        const double& inclinationAngle)
{
    // Declare local variable.
    double pressureCoefficient ;

    // Calculate pressure coefficient.
    pressureCoefficient  = 2 * pow( sin ( inclinationAngle ), 2.0 );
    return pressureCoefficient;
}

//! Function to determine compression and expansion pressure coefficients using
//! Modified Newtonian method.
double calculateModifiedNewtonianPressureCoefficient(
        const double& inclinationAngle,
        const double& stagnationPressureCoefficient)
{
    // Declare local variable.
    double pressureCoefficient;

    // Calculate pressure coefficient.
    pressureCoefficient = stagnationPressureCoefficient *
                          pow ( sin ( inclinationAngle ), 2.0 );
    return pressureCoefficient;
}

//! Function to determine compression pressure coefficients using empirical
//! tangent wedge method.
double  calculateEmpiricalTangentWedgePressureCoefficient(
        const double& inclinationAngle,
        const double& machNumber)
{
    //Declare local variable.
    double machNumberSine;

    //Set local variable.
    machNumberSine = machNumber * sin( inclinationAngle );

    // Form pressure coefficient approximation.
    double pressureCoefficient = ( pow ( 1.2 * machNumberSine +
        exp( -0.6 * machNumberSine ), 2.0 ) - 1 ) /
            ( 0.6 * pow( machNumber, 2.0 ) );

    return pressureCoefficient;
}

//! Function to determine pressure coefficients using empirical
//! tangent cone method.
double  calculateEmpiricalTangentConePressureCoefficient(
        const double& inclinationAngle,
        const double& machNumber)
{
    // Declare local variables.
    double machNumberSine, temp;

    // Set local variables.
    machNumberSine = machNumber * sin ( inclinationAngle );
    temp = pow ( ( 1.090909 * machNumberSine +
                   exp( -0.5454545 * machNumberSine ) ), 2 );

    // Form pressure coefficient approximation.
    double pressureCoefficient = ( 48 * temp *
           pow ( sin ( inclinationAngle), 2.0 ) ) / ( 23 * temp - 5 );

    return pressureCoefficient;
}

//! Function to determine pressure coefficients using modified
//! Dahlem Buck method.
double  calculateModifiedDahlemBuckPressureCoefficient(
        const double& inclinationAngle,
        const double& machNumber)
{
    // Declare local variables.
    double checkAngle = 22.5 * M_PI / 180;
    double factor1;
    double factor2;
    double exponent;
    double pressureCoefficient;

    // Check if inclination angle is greater than check angle. If so, use
    // Newtonian approximation.
    if( inclinationAngle > checkAngle )
    {
        pressureCoefficient = calculateNewtonianPressureCoefficient( inclinationAngle );
    }

    // Else use Dahlem-Buck method.
    else
    {
        pressureCoefficient = ( 1 + sin( 4 * pow( inclinationAngle , 0.75) ) ) /
            ( pow( 4 * cos( inclinationAngle ) * cos( 2 * inclinationAngle ),
                   0.75 ) ) * pow( sin( inclinationAngle ), 1.25 );
    }

    // For mach < 20, a correction term should be applied.
    if( machNumber > 20 )
    {
        factor2 = 1.0;
    }
    else
    {
        // Determine correction term.
        factor1 = ( 6 - 0.3 * machNumber ) + sin( M_PI * (log( machNumber )
                                                          - 0.588 ) / 1.20);
        exponent = 1.15 + 0.5*sin( M_PI * ( log( machNumber ) - 0.916 ) / 3.29 );
        factor2 = 1 + factor1 * pow( inclinationAngle * 180 / M_PI, -1 * exponent );
    }
    // Form pressure coefficient.
    pressureCoefficient = pressureCoefficient * factor2;
    return pressureCoefficient;
}

//! Function to determine compression pressure coefficients using Hankey
//! Flat Surface method.
double  calculateHankeyFlatSurfacePressureCoefficient(
        const double& inclinationAngle,
        const double& machNumber)
{
    // Declare local variables.
    double stagnationPressureCoefficient;
    double pressureCoefficient;

    // Calculate 'effective' stagnation pressure coefficient for low
    // inclination angle.
    if( inclinationAngle < M_PI / 18 )
    {
        stagnationPressureCoefficient = ( 0.195 + 0.222594 /
            pow( machNumber, 0.3 ) - 0.4 ) * inclinationAngle * 180 / M_PI + 4;
    }
    // Calculate 'effective' stagnation pressure coefficient for other
    // inclination angle.
    else
    {
        stagnationPressureCoefficient = 1.95 + 0.3925 /
           ( pow( machNumber, 0.3 ) * tan( inclinationAngle ) );
    }

    // Calculate pressure coefficient using 'effective' stagnation pressure
    // coefficient.
    pressureCoefficient = calculateModifiedNewtonianPressureCoefficient(
            inclinationAngle,
            stagnationPressureCoefficient);
    return pressureCoefficient;
}

//! Function to determine compression pressure coefficients using Smyth
//! Delta Wing method.
double  calculateSmythDeltaWingPressureCoefficient(
        const double& inclinationAngle,
        const double& machNumber)
{
    // Declare local variables.
    double machNumberSine;
    double pressureCoefficient;
    double correctedInclinationAngle;

    // Calculate inclination angle for use in calculations ( angles lower than
    // 1 degree not allowed ).
    if( inclinationAngle<M_PI/180 )
    {
        correctedInclinationAngle = M_PI/180;
    }
    else
    {
        correctedInclinationAngle = inclinationAngle;
    }

    // Pre-compute for efficiency.
    machNumberSine = machNumber * sin( correctedInclinationAngle );

    // Employ empirical correlation to calculate pressure coefficient.
    pressureCoefficient = 1.66667* ( pow (1.09 * machNumberSine +
        exp( -0.49 * machNumberSine ), 2.0 ) - 1 ) / pow( machNumber, 2.0 );

    return pressureCoefficient;
}

//! Function to determine compression and expansion pressure coefficients
//! using Van Dyke Unified method.
double  calculateVanDykeUnifiedPressureCoefficient(
        const double& inclinationAngle,
        const double& machNumber,
        const double& ratioOfSpecificHeats,
        const int& type)
{
    // Declare local variables and pre-compute for efficiency.
    double ratioOfSpecificHeatsTerm = ( ratioOfSpecificHeats + 1 ) / 2;
    double machNumberTerm = pow( pow( machNumber , 2.0 ) - 1, 0.5);
    double exponent = 2 * ratioOfSpecificHeats / ( ratioOfSpecificHeats - 1 );
    double pressureCoefficient;

    // Calculate compression pressure coefficient.
    if( inclinationAngle >= 0 && type == 1 )
    {
        pressureCoefficient = pow( inclinationAngle, 2.0 ) *
            ( ratioOfSpecificHeatsTerm + pow( pow ( ratioOfSpecificHeatsTerm, 2.0 ) +
             4 / ( pow(  inclinationAngle * machNumberTerm, 2.0 ) ), 0.5) );

    }

    //Calculate expansion pressure coefficient.
    else if( inclinationAngle < 0 && type == -1 )
    {
        // Calculate vacuum pressure coefficient.
        double vacuumPressureCoefficient =
                aerodynamics::calculateVacuumPressureCoefficient(
                            machNumber,
                            ratioOfSpecificHeats);

        // Check to see if pressure coefficient will be lower than vacuum case,
        // set to vacuum if so.
        if( -1 * inclinationAngle * machNumberTerm > 2 / ( ratioOfSpecificHeats - 1 ) )
        {
            pressureCoefficient = vacuumPressureCoefficient;
        }
        else
        {
            pressureCoefficient = 2 / ( ratioOfSpecificHeats *
            pow( machNumberTerm, 2.0 ) ) * ( pow( 1 - ( ratioOfSpecificHeats - 1 )
            / 2 * - 1 * inclinationAngle * machNumberTerm, exponent ) - 1 );

            if ( pressureCoefficient < vacuumPressureCoefficient)
            {
                pressureCoefficient = vacuumPressureCoefficient;
            }
        }
    }
    return pressureCoefficient;
}

//! Function to determine expansion pressure coefficients using Prandtl Meyer
//! expansion from freestream.
double  calculatePrandtlMeyerFreestreamPressureCoefficient(
        const double& inclinationAngle,
        const double& machNumber,
        const double& ratioOfSpecificHeats,
        const double& freestreamPrandtlMeyerFunction)
{
    // Declare local variables.
    double prandtlMeyerFunction, pressureCoefficient;

    // Determine Prandtl-Meyer function value.
    prandtlMeyerFunction = freestreamPrandtlMeyerFunction - inclinationAngle;

    // If Prandtl-Meyer function is greater than the vacuum value, set vacuum
    // pressure coefficient.
    if( prandtlMeyerFunction > aerodynamics::maximumPrandtlMeyerFunctionValue )
    {
        pressureCoefficient = aerodynamics::
                              calculateVacuumPressureCoefficient(machNumber,
                                                         ratioOfSpecificHeats);
    }
    else
    {
        // Determine local mach number.
        double localMachNumber =
                aerodynamics::calculateInversePrandtlMeyerFunction(
                                         prandtlMeyerFunction );

        // Determine local to freestream pressure ratio
        double pressureRatio =
                aerodynamics::calculateLocalToStaticPressureRatio(
                        localMachNumber,
                        ratioOfSpecificHeats ) /
                aerodynamics::calculateLocalToStaticPressureRatio(
                        machNumber,
                        ratioOfSpecificHeats );

        // Form pressure coefficient.
        pressureCoefficient = 2 / ( ratioOfSpecificHeats *
                                    pow( machNumber, 2.0 ) ) *
                                    ( pressureRatio -1 ) ;
    }
    return pressureCoefficient;
}

//! Function to calculate pressure coefficient at vacuum.
double  calculateVacuumPressureCoefficient(
        const double& machNumber,
        const double& ratioOfSpecificHeats)
{
    // Declare local variable.
    double pressureCoefficient ;

    // Calculate pressure coefficient.
    pressureCoefficient = -2 / ( ratioOfSpecificHeats * pow( machNumber, 2.0 ) );
    return pressureCoefficient;
}

//! Function to calculate High Mach Base Pressure coefficient.
double  calculateHighMachBasePressure(
        const double& machNumber)
{
    // Declare local variable.
    double pressureCoefficient;

    // Calculate pressure coefficient.
    pressureCoefficient = -1 / pow( machNumber, 2.0 );
    return pressureCoefficient;

}

//! Function to determine expansion pressure coefficients from ACM empirical
//! method.
double  calculateACMempiricalPressureCoefficient(
        const double& inclinationAngle,
        const double& machNumber )
{
    // Declare local variables.
    double pressureCoefficient, minimumPressureCoefficient, testPressureCoefficient;

    // Set minimum pressure coefficient.
    minimumPressureCoefficient = -1 / pow( machNumber, 2.0 );

    // Calculate preliminary pressure coefficient.
    testPressureCoefficient = 180 / M_PI * inclinationAngle / ( 16 * pow(
            machNumber , 2.0 ) );

    // If necessary, correct preliminary pressure coefficient.
    if( minimumPressureCoefficient > testPressureCoefficient )
    {
        pressureCoefficient = minimumPressureCoefficient;
    }
    else
    {
        pressureCoefficient = testPressureCoefficient;
    }

    return pressureCoefficient;
}

//! Function to determine inverse Prandtl-Meyer function for air only.
double  calculateInversePrandtlMeyerFunction(
        const double& prandtlMeyerFunctionValue )
{
    // Declare local variables.
    double y, Mach;

    // Determine input variable for correlation.
    y = pow( prandtlMeyerFunctionValue /
             aerodynamics::maximumPrandtlMeyerFunctionValue, 2.0 / 3.0 );

    // Calculate Mach number.
    Mach = ( 1 + y * ( PrandtlMeyerParameter1 +
                       y * ( PrandtlMeyerParameter2 +
                             y * PrandtlMeyerParameter3 ) ) ) /
            ( 1 + y * ( PrandtlMeyerParameter4 + y * PrandtlMeyerParameter5));

    return Mach;
}

//! Function to calculate ratio of post- to pre-shock pressure.
double calculateShockPressureRatio(
        const double& normalMachNumber,
        const double& ratioOfSpecificHeats)
{
    // Declare local variable.
    double pressureRatio;

    // Calculate pressure ratio.
    pressureRatio = 1 + 2 * ratioOfSpecificHeats / ( ratioOfSpecificHeats + 1 ) *
                                     ( normalMachNumber * normalMachNumber - 1 );

    return pressureRatio;
}

//! Function to calculate ratio of post- to pre-shock density.
double calculateShockDensityRatio(
        const double& normalMachNumber,
        const double& ratioOfSpecificHeats)
{
    // Declare local variables.
    double densityRatio, machNumberSquared ;

    // Calculate mach number squared for efficiency.
    machNumberSquared  = normalMachNumber * normalMachNumber;

    // Calculate density ratio.
    densityRatio = ( ratioOfSpecificHeats + 1 ) * machNumberSquared /
            ( 2 + ( ratioOfSpecificHeats - 1 ) * machNumberSquared );

    return densityRatio;
}

//! Function to calculate ratio of post- to pre-shock temperature.
double calculateShockTemperatureRatio(
        const double& normalMachNumber,
        const double& ratioOfSpecificHeats)
{
    // Declare local variable.
    double temperatureRatio;

    // Calculate temperature ratio from perfect gas law.
    temperatureRatio = 1.0 / aerodynamics::calculateShockDensityRatio(normalMachNumber,
                                                                ratioOfSpecificHeats ) *
                       aerodynamics::calculateShockPressureRatio(normalMachNumber,
                                                                 ratioOfSpecificHeats);

    return temperatureRatio;
}

//! Function to calculate jump in entropy across a shock wave.
double calculateShockEntropyJump(
        const double& normalMachNumber,
        const double& ratioOfSpecificHeats,
        const double& specificGasConstant)
{
    double entropyJump, specificHeatConstantPressure;

    // Calculate specific heat at constant pressure.
    specificHeatConstantPressure = ratioOfSpecificHeats * specificGasConstant /
                                  ( ratioOfSpecificHeats - 1 );

    // Calculate entropy jump from temperature and pressure ratio.
    entropyJump = specificHeatConstantPressure
                  * log( aerodynamics::calculateShockTemperatureRatio(
                          normalMachNumber, ratioOfSpecificHeats ) ) -
                  specificGasConstant *
                  log( aerodynamics::calculateShockPressureRatio(
                          normalMachNumber, ratioOfSpecificHeats ) );
    return entropyJump;
}

//! Function to calculate post- to pre- shock total pressure ratio.
double calculateShockTotalPressureRatio(
            const double& normalMachNumber,
            const double& ratioOfSpecificHeats,
            const double& specificGasConstant)
{
    // Declare local variable.
    double totalPressureRatio;

    // Calculate total pressure ratio from entropy jump.
    totalPressureRatio = exp( -1 * aerodynamics::calculateShockEntropyJump(
            normalMachNumber, ratioOfSpecificHeats, specificGasConstant ) /
                              specificGasConstant);

    return totalPressureRatio;
}

//! Function to calculate shock deflection angle.
double calculateShockDeflectionAngle(
        const double &shockAngle,
        const double &machNumber,
        const double &ratioOfSpecificHeats)
{
    // Declare local variables.
    double tangentOfDeflectionAngle, deflectionAngle;

    // Calculate tangent of deflection angle.
    tangentOfDeflectionAngle = 2 * ( pow( machNumber * sin( shockAngle ) , 2.0 ) - 1 ) /
                               ( tan( shockAngle ) * ( pow( machNumber, 2.0) * (
                               ratioOfSpecificHeats + cos(2*shockAngle)) + 2 ) );

    // Calculate deflection angle.
    deflectionAngle = atan( tangentOfDeflectionAngle );

    return deflectionAngle;
}

}

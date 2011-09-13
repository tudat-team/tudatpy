/*! \file aerodynamics.cpp
 *    This file contains the definition of the aerodynamics namespace.
 *
 *    Path              : /Astrodynamics/ForceModels/Aerothermodynamics/
 *    Version           : 6
 *    Check status      : Checked
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Checker           : L. Abdulkadir
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : L.Abdulkadir@student.tudelft.nl
 *
 *    Date created      : 25 November, 2010
 *    Last modified     : 10 February, 2011
 *
 *    References
 *      Anderson Jr., J.D. , Fundamentals of Aerodynamics, 3rd edition,
 *          McGraw Hill, 2001.
 *      Gentry, A., Smyth, D., and Oliver, W. . The Mark IV Supersonic-
 *          Hypersonic Arbitrary Body Program, Volume II - Program Formulation,
 *          Douglas Aircraft Company, 1973.
 *      Anderson Jr. , J.D, Hypersonic and High-Temperature Gas Dynamics, 2nd
 *          edition, AIAA Education Series, 2006.
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
 *      102511    D. Dirkx          First version of file.
 *      110501    D. Dirkx          Added more comments.
 *      110203    L. Abdulkadir     Code check.
 *      110208    D. Dirkx          Fixed shock temperature ratio.
 *      110210    L. Abdulkadir     Code check.
 *      110211    K. Kumar          Corrected Doxygen errors; corrected layout
 *                                  errors; corrected function-naming;
 *                                  optimized code; corrected double precision.
 */

// Include statements.
#include "aerodynamics.h"

// Using declarations.
using mathematics::raiseToIntegerPower;

namespace aerodynamics
{

//! Compute local-to-static pressure ratio.
double computeLocalToStaticPressureRatio( const double& machNumber,
                                          const double& ratioOfSpecificHeats )
{
    // Return local-to-static pressure ratio.
    return pow( 2.0 / ( 2.0 + ( ratioOfSpecificHeats - 1.0 )
                        * raiseToIntegerPower( machNumber, 2 ) ),
                ratioOfSpecificHeats / ( ratioOfSpecificHeats - 1.0 ) );
}

//! Compute Prandtl-Meyer function.
double computePrandtlMeyerFunction( const double& machNumber,
                                    const double& ratioOfSpecificHeats )
{
    // Declare local variables.
    // Declare Mach number squared.
    double machNumberSquared_ = raiseToIntegerPower( machNumber, 2 );

    // Return value of Prandtl-Meyer function.
    return sqrt ( ( ratioOfSpecificHeats + 1.0 ) /
                  ( ratioOfSpecificHeats - 1.0 ) )
            * atan ( sqrt ( ( ratioOfSpecificHeats - 1.0 ) /
                            ( ratioOfSpecificHeats + 1.0 )
                            * ( machNumberSquared_ - 1.0 ) ) )
            - atan( sqrt ( machNumberSquared_ - 1.0 ) );
}

//! Compute stagnation pressure coefficient in supersonic flow.
double computeStagnationPressure( const double& machNumber,
                                  const double& ratioOfSpecificHeats )
{
    // Declare local variables.
    // Declare Mach number squared.
    double machNumberSquared_ = raiseToIntegerPower( machNumber, 2 );

    // Return stagnation pressure coefficient.
    return 2.0 / ( ratioOfSpecificHeats * machNumberSquared_ )
            * ( pow ( raiseToIntegerPower( ( ratioOfSpecificHeats + 1.0 )
                                           * machNumber, 2 )
                      / ( 4.0 * ratioOfSpecificHeats * machNumberSquared_
                          - 2.0 * ( ratioOfSpecificHeats - 1.0 ) ),
                      ratioOfSpecificHeats / ( ratioOfSpecificHeats - 1.0 ) )
                * ( ( 1.0 - ratioOfSpecificHeats
                      + 2.0 * ratioOfSpecificHeats * machNumberSquared_ )
                    / ( ratioOfSpecificHeats + 1.0 ) ) - 1.0 );
}

//! Compute pressure coefficient based on Newtonian theory.
double computeNewtonianPressureCoefficient( const double& inclinationAngle )
{
    // Return pressure coefficient.
    return 2.0 * raiseToIntegerPower( sin( inclinationAngle ), 2 );
}

//! Compute pressure coefficient based on modified Newtonian theory.
double computeModifiedNewtonianPressureCoefficient(
        const double& inclinationAngle,
        const double& stagnationPressureCoefficient )
{
    // Return pressure coefficient.
    return stagnationPressureCoefficient
            * raiseToIntegerPower( sin( inclinationAngle ), 2 );
}

//! Compute pressure coefficient using empirical tangent wedge method.
double computeEmpiricalTangentWedgePressureCoefficient(
        const double& inclinationAngle, const double& machNumber )
{
    // Declare local variable.
    double machNumberSine_;

    // Set local variable.
    machNumberSine_ = machNumber * sin( inclinationAngle );

    // Return pressure coefficient approximation.
    return ( raiseToIntegerPower( 1.2 * machNumberSine_
                                  + exp( -0.6 * machNumberSine_ ), 2 )
             - 1.0 ) / ( 0.6 * raiseToIntegerPower( machNumber, 2 ) );
}

//! Compute pressure coefficient using empirical tangent cone method.
double computeEmpiricalTangentConePressureCoefficient(
        const double& inclinationAngle, const double& machNumber )
{
    // Declare local variables.
    double machNumberSine_;
    double temporaryValue_;

    // Set local variables.
    machNumberSine_ = machNumber * sin ( inclinationAngle );
    temporaryValue_ = raiseToIntegerPower( ( 1.090909 * machNumberSine_
                                             +  exp( -0.5454545
                                                     * machNumberSine_ ) ), 2 );

    // Return pressure coefficient approximation.
    return ( 48.0 * temporaryValue_
             * raiseToIntegerPower( sin( inclinationAngle), 2 ) )
            / ( 23.0 * temporaryValue_ - 5.0 );
}

//! Compute pressure coefficient using modified Dahlem-Buck method.
double computeModifiedDahlemBuckPressureCoefficient(
        const double& inclinationAngle, const double& machNumber )
{
    // Declare local variables.
    double checkAngle_ = 22.5 * M_PI / 180.0;
    double factor1_;
    double factor2_;
    double exponent_;
    double pressureCoefficient_;

    // Check if inclination angle is greater than check angle. If so, use
    // Newtonian approximation.
    if ( inclinationAngle > checkAngle_ )
    {
        pressureCoefficient_
                = computeNewtonianPressureCoefficient( inclinationAngle );
    }

    // Else use Dahlem-Buck method.
    else
    {
        pressureCoefficient_
                = ( 1.0 + sin( 4.0 * pow( inclinationAngle , 0.75 ) ) )
                  / ( pow( 4.0 * cos( inclinationAngle )
                           * cos( 2.0 * inclinationAngle ), 0.75 ) )
                  * pow( sin( inclinationAngle ), 1.25 );
    }

    // For mach < 20, a correction term should be applied.
    if ( machNumber > 20.0 )
    {
        factor2_ = 1.0;
    }
    else
    {
        // Determine correction term.
        factor1_ = ( 6.0 - 0.3 * machNumber )
                   + sin( M_PI * ( log( machNumber ) - 0.588 ) / 1.20 );

        exponent_ = 1.15
                    + 0.5 * sin( M_PI * ( log( machNumber ) - 0.916 ) / 3.29 );

        factor2_ = 1.0 + factor1_
                   * pow( inclinationAngle * 180.0 / M_PI, -1.0 * exponent_ );
    }

    // Return pressure coefficient.
    return pressureCoefficient_ * factor2_;
}

//! Compute pressure coefficient using the Hankey flat surface method.
double computeHankeyFlatSurfacePressureCoefficient(
        const double& inclinationAngle, const double& machNumber )
{
    // Declare local variables.
    double stagnationPressureCoefficient_;

    // Calculate 'effective' stagnation pressure coefficient for low
    // inclination angle.
    if( inclinationAngle < M_PI / 18.0 )
    {
        stagnationPressureCoefficient_
                = ( 0.195 + 0.222594 / pow( machNumber, 0.3 ) - 0.4 )
                  * inclinationAngle * 180.0 / M_PI + 4.0;
    }
    // Calculate 'effective' stagnation pressure coefficient for other
    // inclination angle.
    else
    {
        stagnationPressureCoefficient_
                = 1.95 + 0.3925 / ( pow( machNumber, 0.3 )
                                    * tan( inclinationAngle ) );
    }

    // Return pressure coefficient using 'effective' stagnation pressure
    // coefficient.
    return computeModifiedNewtonianPressureCoefficient(
            inclinationAngle, stagnationPressureCoefficient_ );
}

//! Compute pressure coefficient using the Smyth delta wing method.
double computeSmythDeltaWingPressureCoefficient(
        const double& inclinationAngle, const double& machNumber )
{
    // Declare local variables.
    double machNumberSine_;
    double correctedInclinationAngle_;

    // Calculate inclination angle for use in calculations ( angles lower than
    // 1 degree not allowed ).
    if ( inclinationAngle < M_PI / 180 )
    {
        correctedInclinationAngle_ = M_PI / 180;
    }

    else
    {
        correctedInclinationAngle_ = inclinationAngle;
    }

    // Pre-compute for efficiency.
    machNumberSine_ = machNumber * sin( correctedInclinationAngle_ );

    // Employ empirical correlation to calculate pressure coefficient.
    // Return pressure coefficient.
    return 1.66667 * ( raiseToIntegerPower(
            1.09 * machNumberSine_ + exp( -0.49 * machNumberSine_ ), 2 ) - 1.0 )
            / raiseToIntegerPower( machNumber, 2 );
}

//! Compute pressure coefficient using the van Dyke unified method.
double computeVanDykeUnifiedPressureCoefficient(
        const double& inclinationAngle, const double& machNumber,
        const double& ratioOfSpecificHeats, const int& type )
{
    // Declare and initialize local variables and pre-compute for efficiency.
    double ratioOfSpecificHeatsTerm_ = ( ratioOfSpecificHeats + 1.0 ) / 2.0;
    double machNumberTerm_ = sqrt( raiseToIntegerPower( machNumber , 2 )
                                   - 1.0 );
    double exponent_ = 2.0
                       * ratioOfSpecificHeats / ( ratioOfSpecificHeats - 1.0 );

    // Declare and initialize value.
    double pressureCoefficient_ = 0.0;

    // Calculate compression pressure coefficient.
    if ( inclinationAngle >= 0 && type == 1 )
    {
        pressureCoefficient_
                = raiseToIntegerPower( inclinationAngle, 2 )
                  * ( ratioOfSpecificHeatsTerm_
                      + sqrt( raiseToIntegerPower( ratioOfSpecificHeatsTerm_, 2 )
                              + 4.0
                              / ( raiseToIntegerPower(  inclinationAngle
                                                        * machNumberTerm_,
                                                        2 ) ) ) );
    }

    // Calculate expansion pressure coefficient.
    else if ( inclinationAngle < 0 && type == -1 )
    {
        // Calculate vacuum pressure coefficient.
        double vacuumPressureCoefficient_ = computeVacuumPressureCoefficient(
                machNumber, ratioOfSpecificHeats );

        // Check to see if pressure coefficient will be lower than vacuum case,
        // set to vacuum if so.
        if ( -1.0 * inclinationAngle * machNumberTerm_
             > 2.0 / ( ratioOfSpecificHeats - 1.0 ) )
        {
            pressureCoefficient_ = vacuumPressureCoefficient_;
        }
        else
        {
            pressureCoefficient_
                    = 2.0 / ( ratioOfSpecificHeats
                              * raiseToIntegerPower( machNumberTerm_, 2 ) )
                    * ( pow( 1.0 - ( ratioOfSpecificHeats - 1.0 ) / 2.0
                             * - 1.0 * inclinationAngle
                             * machNumberTerm_, exponent_ ) - 1.0 );

            if ( pressureCoefficient_ < vacuumPressureCoefficient_ )
            {
                pressureCoefficient_ = vacuumPressureCoefficient_;
            }
        }
    }

    // Return pressure coefficient.
    return pressureCoefficient_;
}

//! Compute pressure coefficient using Prandtl-Meyer expansion.
double computePrandtlMeyerFreestreamPressureCoefficient(
        const double& inclinationAngle,  const double& machNumber,
        const double& ratioOfSpecificHeats,
        const double& freestreamPrandtlMeyerFunction )
{
    // Declare local variables.
    double prandtlMeyerFunction_;
    double pressureCoefficient_;

    // Determine Prandtl-Meyer function value.
    prandtlMeyerFunction_ = freestreamPrandtlMeyerFunction - inclinationAngle;

    // If Prandtl-Meyer function is greater than the vacuum value, set vacuum
    // pressure coefficient.
    if ( prandtlMeyerFunction_ > maximumPrandtlMeyerFunctionValue )
    {
        pressureCoefficient_ = computeVacuumPressureCoefficient(
                machNumber, ratioOfSpecificHeats );
    }

    else
    {
        // Determine local mach number.
        double localMachNumber_
                = computeInversePrandtlMeyerFunction( prandtlMeyerFunction_ );

        // Determine local to freestream pressure ratio.
        double pressureRatio_
                = computeLocalToStaticPressureRatio( localMachNumber_,
                                                     ratioOfSpecificHeats )
                / computeLocalToStaticPressureRatio( machNumber,
                                                     ratioOfSpecificHeats );

        // Form pressure coefficient.
        pressureCoefficient_ = 2.0 / ( ratioOfSpecificHeats
                                       *  raiseToIntegerPower( machNumber, 2 ) )
                               * ( pressureRatio_ - 1.0 ) ;
    }

    // Return pressure coefficient.
    return pressureCoefficient_;
}

//! Compute pressure coefficient at vacuum.
double computeVacuumPressureCoefficient(
        const double& machNumber, const double& ratioOfSpecificHeats )
{
    // Return pressure coefficient.
    return -2.0 / ( ratioOfSpecificHeats
                    * raiseToIntegerPower( machNumber, 2 ) );
}

//! Compute high Mach base pressure coefficient.
double computeHighMachBasePressure( const double& machNumber )
{
    // Calculate pressure coefficient.
    return -1.0 / raiseToIntegerPower( machNumber, 2 );
}

//! Compute pressure coefficient using the ACM empirical method.
double computeAcmEmpiricalPressureCoefficient(
        const double& inclinationAngle, const double& machNumber )
{
    // Declare local variables.
    double pressureCoefficient_;
    double minimumPressureCoefficient_;
    double preliminaryPressureCoefficient_;

    // Set minimum pressure coefficient.
    minimumPressureCoefficient_ = -1.0 / raiseToIntegerPower( machNumber, 2 );

    // Calculate preliminary pressure coefficient.
    preliminaryPressureCoefficient_ = 180.0 / M_PI * inclinationAngle
                                      / ( 16.0 * raiseToIntegerPower(
                                              machNumber , 2 ) );

    // If necessary, correct preliminary pressure coefficient.
    if ( minimumPressureCoefficient_ > preliminaryPressureCoefficient_ )
    {
        pressureCoefficient_ = minimumPressureCoefficient_;
    }

    else
    {
        pressureCoefficient_ = preliminaryPressureCoefficient_;
    }

    // Return pressure coefficient.
    return pressureCoefficient_;
}

//! Compute Mach number from Prandtl-Meyer function.
double computeInversePrandtlMeyerFunction( const double&
                                           prandtlMeyerFunctionValue )
{
    // Declare local variables.
    double inputVariableForCorrelation_;
    double machNumber_;

    // Determine input variable for correlation.
    inputVariableForCorrelation_ = pow( prandtlMeyerFunctionValue
                                        / maximumPrandtlMeyerFunctionValue,
                                        2.0 / 3.0 );

    // Calculate Mach number.
    machNumber_ = ( 1.0 + inputVariableForCorrelation_
                    * ( PrandtlMeyerParameter1 + inputVariableForCorrelation_
                        * ( PrandtlMeyerParameter2
                            + inputVariableForCorrelation_
                            * PrandtlMeyerParameter3 ) ) )
            / ( 1.0 + inputVariableForCorrelation_
                * ( PrandtlMeyerParameter4
                    + inputVariableForCorrelation_
                    * PrandtlMeyerParameter5 ) );

    // Return Mach number.
    return machNumber_;
}

//! Compute ratio of post- to pre-shock pressure.
double computeShockPressureRatio( const double& normalMachNumber,
                                  const double& ratioOfSpecificHeats )
{
    // Return pressure ratio.
    return 1.0 + 2.0 * ratioOfSpecificHeats
            / ( ratioOfSpecificHeats + 1.0 )
            * ( normalMachNumber * normalMachNumber - 1.0 );
}

//! Compute ratio of post- to pre-shock density.
double computeShockDensityRatio( const double& normalMachNumber,
                                 const double& ratioOfSpecificHeats )
{
    // Declare local variables.
    double machNumberSquared_;

    // Calculate mach number squared for efficiency.
    machNumberSquared_ = raiseToIntegerPower( normalMachNumber, 2 );

    // Return density ratio.
    return ( ratioOfSpecificHeats + 1.0 ) * machNumberSquared_
            / ( 2.0 + ( ratioOfSpecificHeats - 1.0 ) * machNumberSquared_ );
}

//! Compute ratio of post- to pre-shock temperature.
double computeShockTemperatureRatio( const double& normalMachNumber,
                                     const double& ratioOfSpecificHeats )
{
    // Return temperature ratio from perfect gas law.
    return 1.0 / computeShockDensityRatio( normalMachNumber,
                                           ratioOfSpecificHeats )
           * computeShockPressureRatio( normalMachNumber,
                                        ratioOfSpecificHeats );
}

//! Compute jump in entropy across a shock wave.
double computeShockEntropyJump( const double& normalMachNumber,
                                const double& ratioOfSpecificHeats,
                                const double& specificGasConstant )
{
    // Declare local variables.
    double specificHeatConstantPressure_;

    // Calculate specific heat at constant pressure.
    specificHeatConstantPressure_ = ratioOfSpecificHeats
                                    * specificGasConstant
                                    / ( ratioOfSpecificHeats - 1.0 );

    // Return entropy jump from temperature and pressure ratio.
    return specificHeatConstantPressure_
            * log( computeShockTemperatureRatio( normalMachNumber,
                                                 ratioOfSpecificHeats ) )
            - specificGasConstant
            * log( computeShockPressureRatio( normalMachNumber,
                                              ratioOfSpecificHeats ) );
}

//! Compute post- to pre-shock total pressure ratio.
double computeShockTotalPressureRatio( const double& normalMachNumber,
                                       const double& ratioOfSpecificHeats,
                                       const double& specificGasConstant )
{
    // Return total pressure ratio from entropy jump.
    return exp( -1.0 * computeShockEntropyJump( normalMachNumber,
                                                ratioOfSpecificHeats,
                                                specificGasConstant )
                / specificGasConstant );
}

//! Compute shock deflection angle.
double computeShockDeflectionAngle( const double& shockAngle,
                                    const double& machNumber,
                                    const double& ratioOfSpecificHeats )
{
    // Declare local variables.
    double tangentOfDeflectionAngle_;

    // Calculate tangent of deflection angle.
    tangentOfDeflectionAngle_ = 2.0
                                * ( raiseToIntegerPower( machNumber
                                                         * sin( shockAngle ),
                                                         2 ) - 1.0 )
                                / ( tan( shockAngle )
                                    * ( raiseToIntegerPower( machNumber, 2 )
                                        * ( ratioOfSpecificHeats
                                            + cos( 2.0 * shockAngle ) )
                                        + 2.0 ) );

    // Return deflection angle.
    return atan( tangentOfDeflectionAngle_ );
}

}

// End of file.

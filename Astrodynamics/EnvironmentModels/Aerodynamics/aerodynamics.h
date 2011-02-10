/*! \file aerodynamics.h
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
 *  Last modified     : 08 February, 2011
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
 */

#ifndef AERODYNAMICS_H
#define AERODYNAMICS_H

#include <cmath>

namespace aerodynamics{

    //! Function to calculate local to static pressure ratio.
    /*!
     *  This function calculates the local to static pressure ratio,
     *  assuming a thermally and calorically perfect gas.
     *  \param machNumber Flow mach number.
     *  \param ratioOfSpecificHeats Ratio of specific heat at constant
     *  pressure to specific heat at constant volume.
     *  \return local to static pressure ratio.
     */
    double calculateLocalToStaticPressureRatio(
            const double& machNumber,
            const double& ratioOfSpecificHeats);

    //! Function to evaluate the Prandtl-Meyer function.
    /*!
     *  This functionn calculates the value of the Prandtl-Meyer function
     *  at the given mach number and ratio of specific heat.
     *  \param machNumber Flow mach number.
     *  \param ratioOfSpecificHeats Ratio of specific heat at constant
     *  pressure to specific heat at constant volume.
     *  \return Prandtl-Meyer function value.
     */
    double evaluatePrandtlMeyerFunction(
            const double& machNumber,
            const double& ratioOfSpecificHeats );

    //! Function to calculate stagnation pressure coefficient in supersonic flow.
    /*!
     *  This function calculates the stagnation pressure coefficient,
     *  assuming a thermally and calorically perfect gas.
     *  \param machNumber Flow mach number.
     *  \param ratioOfSpecificHeats Ratio of specific heat at constant
     *  pressure to specific heat at constant volume.
     *  \return stagnation pressure coefficient.
     */
    double calculateStagnationPressure(
            const double& machNumber,
            const double& ratioOfSpecificHeats);

    //! Function to calculate pressure coefficient based on Newtonian theory.
    /*!
     *  This function calculates the pressure coefficient based on Newtonian
     *  theory.
     *  \param inclinationAngle Angle between the wall and the freestream
     *  velocity vector.
     *  \return Newtonian pressure coefficient.
     */
    double calculateNewtonianPressureCoefficient(
            const double& inclinationAngle);

    //! Function to calculate pressure coefficient based on modified Newtonian
    //! theory.
    /*!
     *  This function calculates the pressure coefficient based on modified
     *  Newtonian theory.
     *  \param inclinationAngle Angle between the wall and the freestream
     *  velocity vector.
     *  \param machNumber Flow mach number.
     *  \param ratioOfSpecificHeats ratio of specific heat at constant
     *  pressure to specific heat at constant volume.
     *  \param stagnationPressureCoefficient stagnation pressure coefficient.
     *  \return Newtonian pressure coefficient.
     */
    double calculateModifiedNewtonianPressureCoefficient(
            const double& inclinationAngle,
            const double& stagnationPressureCoefficient);

    //! Function to determine pressure coefficients using empirical tangent
    //! wedge method.
    /*! Function to determine tangent wedge pressure coefficient based on
     *  empirical correlation for ratio of specific heats = 1.4 ( terrestrial
     *  atmosphere).
     *  \param inclinationAngle Angle between the wall and the freestream
     *  velocity vector.
     *  \param machNumber Flow mach number.
     *  \return Empirical tangent wedge pressure coefficient.
     */
    double calculateEmpiricalTangentWedgePressureCoefficient(
            const double& inclinationAngle,
            const double& machNumber);

    //! Function to determine pressure coefficients using empirical tangent
    //! cone method.
    /*! Function to determine tangent cone pressure coefficient based on
     *  empirical correlation for ratio of specific heats = 1.4 ( terrestrial
     *  atmosphere).
     *  \param inclinationAngle Angle between the wall and the freestream
     *  velocity vector.
     *  \param machNumber Flow mach number.
     *  \return Empirical tangent wedge pressure coefficient.
     */
    double calculateEmpiricalTangentConePressureCoefficient(
            const double& inclinationAngle,
            const double& machNumber);

    //! Function to determine pressure coefficients using modified
    //! Dahlem Buck method.
    /*! Function to determine tangent cone pressure coefficient based on
     *  Dahlem-Buck empirical method.
     *  \param inclinationAngle Angle between the wall and the freestream
     *  velocity vector.
     *  \param machNumber Flow mach number.
     *  \return Dahlem-Buck pressure coefficient.
     */
    double calculateModifiedDahlemBuckPressureCoefficient(
            const double& inclinationAngle,
            const double& machNumber);

    //! Function to determine pressure coefficients using the Hankey flat
    //! surface method.
    /*! Function to determine tangent cone pressure coefficient based on
     *  the Hankey flat surface method.
     *  \param inclinationAngle Angle between the wall and the freestream
     *  velocity vector.
     *  \param machNumber Flow mach number.
     *  \return Hankey Flat surface pressure coefficient.
     */
    double calculateHankeyFlatSurfacePressureCoefficient(
            const double& inclinationAngle,
            const double& machNumber);

    //! Function to determine pressure coefficients using the Smyth delta
    //! wing method.
    /*! Function to determine tangent cone pressure coefficient based on
     *  the Smyth delta wing surface method.
     *  \param inclinationAngle Angle between the wall and the freestream
     *  velocity vector.
     *  \param machNumber Flow mach number.
     *  \return Smyth delta wing pressure coefficient.
     */
    double calculateSmythDeltaWingPressureCoefficient(
            const double& inclinationAngle,
            const double& machNumber);

    //! Function to determine pressure coefficients using the van Dyke unified
    //! method.
    /*! Function to determine tangent cone pressure coefficient based on
     *  the van Dyke unified method.
     *  \param inclinationAngle Angle between the wall and the freestream
     *  velocity vector.
     *  \param machNumber Flow mach number.
     *  \param ratioOfSpecificHeats ratio of specific heat at constant
     *  pressure to specific heat at constant volume.
     *  \param type Expansion (1) or compression(-1)
     *  \return Hankey Flat surface pressure coefficient.
     */
    double calculateVanDykeUnifiedPressureCoefficient(
            const double& inclinationAngle,
            const double& machNumber,
            const double& ratioOfSpecificHeats,
            const int& type);

    //! Function to determine pressure coefficients using Prandtl-Meyer
    //! expansion.
    /*! Function to determine pressure coefficients using Prandtl-Meyer
     *  expansion from freestream. Currently only terrestrial atmosphere
     *  ( ratio of specific heat = 1.4 ) is supported due to the use of an
     *  empirical fit for the inverse Prandtl-Meyer function determination.
     *  \param inclinationAngle Angle between the wall and the freestream
     *  velocity vector.
     *  \param machNumber Flow mach number.
     *  \param ratioOfSpecificHeats ratio of specific heat at constant
     *  pressure to specific heat at constant volume.
     *  \param freestreamPrandtlMeyerFunction Freestream Prandtl-Meyer function.
     *  \return Prandtl - Meyer pressure coefficient.
     */
    double calculatePrandtlMeyerFreestreamPressureCoefficient(
            const double& inclinationAngle,
            const double& machNumber,
            const double& ratioOfSpecificHeats,
            const double& freestreamPrandtlMeyerFunction);

    //! Function to calculate pressure coefficient at vacuum.
    /*!
     *  This function calculates the pressure coefficient at vacuum assuming
     *  a thermally and calorically perfect gas.
     *  \param machNumber Flow mach number.
     *  \param ratioOfSpecificHeats ratio of specific heat at constant
     *  pressure to specific heat at constant volume.
     *  \return Vacuum pressure coefficient.
     */
    double calculateVacuumPressureCoefficient(
            const double& machNumber,
            const double& ratioOfSpecificHeats);

    //! Function to calculate High Mach Base Pressure coefficient.
    /*!
     *  This function calculates the high mach base pressure coefficientapproximation
     *  \param machNumber Flow mach number.
     *  \return Vacuum pressure coefficient.
     */
    double calculateHighMachBasePressure(
            const double& machNumber);

    //! Function to determine pressure coefficients using the ACM empirical
    //! method.
    /*! Function to determine tangent cone pressure coefficient based on
     *  the ACM empirical method.
     *  \param inclinationAngle Angle between the wall and the freestream
     *  velocity vector.
     *  \param machNumber Flow mach number.
     *  \return ACM empirical surface pressure coefficient.
     */
    double calculateACMempiricalPressureCoefficient(
            const double& inclinationAngle,
            const double& machNumber );

    //! Function to calculate mach number from Prandtl-Meter function.
    /*!
     *  This function calculates the inverse of the Prandtl-Meyer function.
     *  Currently the function is limited to use with ratio of specific heats
     *  = 1.4 as it uses an empirical correlation.
     *  \param prandtlMeyerFunctionValue Prandyl-Meyer function value.
     *  \return mach number.
     */
    double calculateInversePrandtlMeyerFunction(
            const double& prandtlMeyerFunctionValue);

    //! Function to calculate ratio of post- to pre-shock pressure.
    /*! Function to calculate ratio of post- to pre-shock pressure, assuming
     *  thermally and calorically perfect gas.
     *  \param normalMachNumber Mach number of flow velociry normal to shock.
     *  \param ratioOfSpecificHeats ratio of specific heat at constant
     *  pressure to specific heat at constant volume.
     */
    double calculateShockPressureRatio(
            const double& normalMachNumber,
            const double& ratioOfSpecificHeats);

    //! Function to calculate ratio of post- to pre-shock density.
    /*! Function to calculate ratio of post- to pre-shock density, assuming
     *  thermally and calorically perfect gas.
     *  \param normalMachNumber Mach number of flow velociry normal to shock.
     *  \param ratioOfSpecificHeats ratio of specific heat at constant
     *  pressure to specific heat at constant volume.
     */
    double calculateShockDensityRatio(
            const double& normalMachNumber,
            const double& ratioOfSpecificHeats);


    //! Function to calculate ratio of post- to pre-shock temperature.
    /*! Function to calculate ratio of post- to pre-shock temperature, assuming
     *  thermally and calorically perfect gas.
     *  \param normalMachNumber Mach number of flow velociry normal to shock.
     *  \param ratioOfSpecificHeats ratio of specific heat at constant
     *  pressure to specific heat at constant volume.
     */
    double calculateShockTemperatureRatio(
            const double& normalMachNumber,
            const double& ratioOfSpecificHeats);

    //! Function to calculate jump in entropy across a shock wave.
    /*! Function to calculate jump in entropy across a shock wave, assuming
     *  thermally and calorically perfect gas.
     *  \param normalMachNumber Mach number of flow velociry normal to shock.
     *  \param ratioOfSpecificHeats ratio of specific heat at constant
     *  pressure to specific heat at constant volume.
     *  \param specificGasConstant gas constant per unit mass for flow
     *  composition.
     */
    double calculateShockEntropyJump(
            const double& normalMachNumber,
            const double& ratioOfSpecificHeats,
            const double& specificGasConstant);

    //! Function to calculate post- to pre-shock total pressure ratio.
    /*!
     *  Function to calculate post- to pre-shock total pressure ratio from
     *  the entropy jump across a shock wave. Assumption of thermally and
     *  calorically perfect gas is made.
     *  \param normalMachNumber Mach number of flow velociry normal to shock.
     *  \param ratioOfSpecificHeats ratio of specific heat at constant
     *  pressure to specific heat at constant volume.
     *  \param specificGasConstant gas constant per unit mass for flow
     *  composition.
     */
    double calculateShockTotalPressureRatio(
            const double& normalMachNumber,
            const double& ratioOfSpecificHeats,
            const double& specificGasConstant);

    //! Function to calculate shock deflection angle.
    /*!
     *  Function to calculate the flow deflection angle across a shock wave.
     *  \param shockAngle Angle of shock wave w.r.t. freestream flow.
     *  \param machNumber Freestream Mach number.
     *  \param ratioOfSpecificHeats ratio of specific heat at constant
     *  pressure to specific heat at constant volume.
     */
    double calculateShockDeflectionAngle(
            const double &shockAngle,
            const double &machNumber,
            const double &ratioOfSpecificHeats);

    //! Maximum Prandtl-Meyer function value.
    /*!
     *  Maximum Prandtl-Meyer function value for ratio of specific heats = 1.4.
     */
    static const double maximumPrandtlMeyerFunctionValue = M_PI / 2 *
                                                ( pow( 6, 0.5 ) - 1 );

    //! Constant for use in inverse Prandtl-Meyer function calculation.
    /*!
     *  Constant for use in inverse Prandtl-Meyer function calculation for
     *  ratio of specific heats = 1.4.
     */
    static const double PrandtlMeyerParameter1 = 1.3604;

    //! Constant for use in inverse Prandtl-Meyer function calculation.
    /*!
     *  Constant for use in inverse Prandtl-Meyer function calculation for
     *  ratio of specific heats = 1.4.
     */
    static const double PrandtlMeyerParameter2 = 0.0962;

    //! Constant for use in inverse Prandtl-Meyer function calculation.
    /*!
     *  Constant for use in inverse Prandtl-Meyer function calculation for
     *  ratio of specific heats = 1.4.
     */
    static const double PrandtlMeyerParameter3 = -0.5127;

    //! Constant for use in inverse Prandtl-Meyer function calculation.
    /*!
     *  Constant for use in inverse Prandtl-Meyer function calculation for
     *  ratio of specific heats = 1.4.
     */
    static const double PrandtlMeyerParameter4 = -0.6722;

    //! Constant for use in inverse Prandtl-Meyer function calculation.
    /*!
     *  Constant for use in inverse Prandtl-Meyer function calculation for
     *  ratio of specific heats = 1.4.
     */
    static const double PrandtlMeyerParameter5 = -0.3278;

}

#endif // AERODYNAMICS_H

// End of file.

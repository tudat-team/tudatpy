/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Anderson Jr., J.D. , Fundamentals of aerodynamics, 3rd edition, McGraw Hill, 2001.
 *      Gentry, A., Smyth, D., and Oliver, W. . The Mark IV Supersonic-Hypersonic Arbitrary Body
 *          Program, Volume II - Program Formulation, Douglas Aircraft Company, 1973.
 *      Anderson Jr. , J.D, Hypersonic and High-Temperature Gas Dynamics, 2nd edition,
 *          AIAA Education Series, 2006.
 *
 */

#ifndef TUDAT_AERODYNAMICS_H
#define TUDAT_AERODYNAMICS_H

#include <functional>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <vector>
#include <cmath>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/basics/basicTypedefs.h"

namespace tudat
{
namespace aerodynamics
{

//! Enum defining a list of independent variables on which the aerodynamic coefficients can depend.
/*!
 *  Enum defining a list of independent variables on which the aerodynamic coefficients can depend.
 *  Note that for a custom coefficient interface with other variables, you may use the
 *  undefined_independent_variable variable type, but at the expense of being able to use the
 *  FlightConditions class to automatically updates the aerodynamic coefficients during propagation.
 */
enum AerodynamicCoefficientsIndependentVariables
{
    mach_number_dependent = 0,
    angle_of_attack_dependent = 1,
    angle_of_sideslip_dependent = 2,
    altitude_dependent = 3,
    time_dependent = 4,
    temperature_dependent = 5,
    velocity_dependent = 6,
    he_number_density_dependent = 7,
    o_number_density_dependent = 8,
    n2_number_density_dependent = 9,
    o2_number_density_dependent = 10,
    ar_number_density_dependent = 11,
    h_number_density_dependent = 12,
    n_number_density_dependent = 13,
    anomalous_o_number_density_dependent = 14,
    control_surface_deflection_dependent = 15,
    undefined_independent_variable = 16
};
 
//! Function to combined the force and moment coefficients from separate function pointers.
/*!
 *  Function to combined the force and moment coefficients from separate function pointers.
 *  The output is the concatenated force and moment coefficient vector, evaluated
 *  at the current set of independent variables.
 *  \param forceCoefficientFunction Function returning the aerodynamic force coefficients as
 *  function of the set of independent variables.
 *  \param momentCoefficientFunction Function returning the aerodynamic force coefficients as
 *  function of the set of independent variables.
 *  \param independentVariables Current list of values of the independent variables upon
 *  which the coefficients depend.
 */
inline Eigen::Vector6d concatenateForceAndMomentCoefficients(
        const std::function< Eigen::Vector3d( const std::vector< double >& ) >&
        forceCoefficientFunction,
        const std::function< Eigen::Vector3d( const std::vector< double >& ) >&
        momentCoefficientFunction,
        const std::vector< double >& independentVariables )
{
    return ( Eigen::Vector6d( ) << forceCoefficientFunction( independentVariables ),
             momentCoefficientFunction( independentVariables ) ).finished( );
}

//! Maximum Prandtl-Meyer function value.
/*!
 * Maximum Prandtl-Meyer function value for ratio of specific heats = 1.4.
 */
static const double maximumPrandtlMeyerFunctionValue = 
    mathematical_constants::PI / 2.0 * ( std::sqrt( 6.0 ) - 1.0 );

//! Constant for use in inverse Prandtl-Meyer function calculation.
/*!
 * Constant for use in inverse Prandtl-Meyer function calculation for ratio of
 *  specific heats = 1.4.
 */
static const double PrandtlMeyerParameter1 = 1.3604;

//! Constant for use in inverse Prandtl-Meyer function calculation.
/*!
 * Constant for use in inverse Prandtl-Meyer function calculation for ratio of
 * specific heats = 1.4.
 */
static const double PrandtlMeyerParameter2 = 0.0962;

//! Constant for use in inverse Prandtl-Meyer function calculation.
/*!
 * Constant for use in inverse Prandtl-Meyer function calculation for ratio of
 * specific heats = 1.4.
 */
static const double PrandtlMeyerParameter3 = -0.5127;

//! Constant for use in inverse Prandtl-Meyer function calculation.
/*!
 * Constant for use in inverse Prandtl-Meyer function calculation for ratio of
 * specific heats = 1.4.
 */
static const double PrandtlMeyerParameter4 = -0.6722;

//! Constant for use in inverse Prandtl-Meyer function calculation.
/*!
 * Constant for use in inverse Prandtl-Meyer function calculation for ratio of
 * specific heats = 1.4.
 */
static const double PrandtlMeyerParameter5 = -0.3278;

//! Compute local-to-static pressure ratio.
/*!
 * Computes the local to static pressure ratio, assuming a thermally and
 * calorically perfect gas.
 * \param machNumber Flow Mach number.
 * \param ratioOfSpecificHeats Ratio of specific heat at constant pressure to
 *         specific heat at constant volume.
 * \return Local-to-static pressure ratio.
 */
double computeLocalToStaticPressureRatio( double machNumber,
                                          double ratioOfSpecificHeats );

//! Compute Prandtl-Meyer function.
/*!
 * Computes the value of the Prandtl-Meyer function at the given Mach number
 * and ratio of specific heat.
 * \param machNumber Flow Mach number.
 * \param ratioOfSpecificHeats Ratio of specific heat at constant pressure to
 *         specific heat at constant volume.
 * \return Prandtl-Meyer function value.
 */
double computePrandtlMeyerFunction( double machNumber, double ratioOfSpecificHeats );

//! Compute stagnation pressure coefficient in supersonic flow.
/*!
 * Computes the stagnation pressure coefficient, assuming a thermally and
 * calorically perfect gas.
 * \param machNumber Flow Mach number.
 * \param ratioOfSpecificHeats Ratio of specific heat at constant pressure to
 *         specific heat at constant volume.
 * \return Stagnation pressure coefficient.
 */
double computeStagnationPressure( double machNumber, double ratioOfSpecificHeats );

//! Compute pressure coefficient based on Newtonian theory.
/*!
 * Computes the pressure coefficient based on Newtonian theory.
 * \param inclinationAngle Angle between wall and freestream velocity vector.
 * \return Newtonian pressure coefficient.
 */
double computeNewtonianPressureCoefficient( double inclinationAngle );

//! Compute pressure coefficient based on modified Newtonian theory.
/*!
 * Computes the pressure coefficient based on modified Newtonian theory.
 * \param inclinationAngle Angle between the wall and the freestream
 *         velocity vector.
 * \param stagnationPressureCoefficient Stagnation pressure coefficient.
 * \return Newtonian pressure coefficient.
 */
double computeModifiedNewtonianPressureCoefficient(
        double inclinationAngle, double stagnationPressureCoefficient );

//! Compute pressure coefficient using empirical tangent wedge method.
/*!
 * Computes tangent wedge pressure coefficient based on empirical correlation
 * for ratio of specific heats = 1.4 ( terrestrial atmosphere).
 * \param inclinationAngle Angle between wall and freestream velocity vector.
 * \param machNumber Flow Mach number.
 * \return Empirical tangent wedge pressure coefficient.
 */
double computeEmpiricalTangentWedgePressureCoefficient(
        double inclinationAngle, double machNumber );

//! Compute pressure coefficient using empirical tangent cone method.
/*!
 * Computes tangent cone pressure coefficient based on empirical correlation
 * for ratio of specific heats = 1.4 ( terrestrial atmosphere).
 * \param inclinationAngle Angle between wall and freestream velocity vector.
 * \param machNumber Flow Mach number.
 * \return Empirical tangent wedge pressure coefficient.
 */
double computeEmpiricalTangentConePressureCoefficient(
        double inclinationAngle, double machNumber );

//! Compute pressure coefficient using modified Dahlem-Buck method.
/*!
 * Computes tangent cone pressure coefficient based on Dahlem-Buck empirical
 * method.
 * \param inclinationAngle Angle between wall and freestream velocity vector.
 * \param machNumber Flow Mach number.
 * \return Dahlem-Buck pressure coefficient.
 */
double computeModifiedDahlemBuckPressureCoefficient(
        double inclinationAngle, double machNumber );

//! Compute pressure coefficient using the Hankey flat surface method.
/*!
 * Computes tangent cone pressure coefficient based on the Hankey flat surface
 * method.
 * \param inclinationAngle Angle between wall and freestream velocity vector.
 * \param machNumber Flow Mach number.
 * \return Hankey Flat surface pressure coefficient.
 */
double computeHankeyFlatSurfacePressureCoefficient(
        double inclinationAngle, double machNumber );

//! Compute pressure coefficient using the Smyth delta wing method.
/*!
 * Computes tangent cone pressure coefficient based on the Smyth delta wing
 * surface method.
 * \param inclinationAngle Angle between wall and freestream velocity vector.
 * \param machNumber Flow Mach number.
 * \return Smyth delta wing pressure coefficient.
 */
double computeSmythDeltaWingPressureCoefficient(
        double inclinationAngle, double machNumber );

//! Compute pressure coefficient using the van Dyke unified method.
/*!
 * Computes tangent cone pressure coefficient based on the van Dyke unified
 * method.
 * \param inclinationAngle Angle between wall and freestream velocity vector.
 * \param machNumber Flow Mach number.
 * \param ratioOfSpecificHeats Ratio of specific heat at constant pressure to
 *         specific heat at constant volume.
 * \param type ( expansion ( 1 ) or compression( -1 ) ).
 * \return Hankey Flat surface pressure coefficient.
 */
double computeVanDykeUnifiedPressureCoefficient(
        double inclinationAngle, double machNumber,
        double ratioOfSpecificHeats, int type );

//! Compute pressure coefficient using Prandtl-Meyer expansion.
/*!
 * Computes pressure coefficient using Prandtl-Meyer expansion from
 * freestream. Currently only terrestrial atmosphere
 * ( ratio of specific heat = 1.4 ) is supported due to the use of an empirical
 * fit for the inverse Prandtl-Meyer function determination.
 * \param inclinationAngle Angle between wall and freestream velocity vector.
 * \param machNumber Flow Mach number.
 * \param ratioOfSpecificHeats Ratio of specific heat at constant pressure to
 *         specific heat at constant volume.
 * \param freestreamPrandtlMeyerFunction Freestream Prandtl-Meyer function.
 * \return Prandtl-Meyer pressure coefficient.
 */
double computePrandtlMeyerFreestreamPressureCoefficient(
        double inclinationAngle,  double machNumber,
        double ratioOfSpecificHeats, double freestreamPrandtlMeyerFunction );

//! Compute pressure coefficient at vacuum.
/*!
 * Computes the pressure coefficient at vacuum assuming a thermally and
 * calorically perfect gas.
 * \param machNumber Flow Mach number.
 * \param ratioOfSpecificHeats Ratio of specific heat at constant pressure to
 *         specific heat at constant volume.
 * \return Vacuum pressure coefficient.
 */
double computeVacuumPressureCoefficient(
        double machNumber, double ratioOfSpecificHeats );

//! Compute high Mach base pressure coefficient.
/*!
 * This function calculates the high Mach base pressure coefficient
 * approximation.
 * \param machNumber Flow Mach number.
 * \return Vacuum pressure coefficient.
 */
double computeHighMachBasePressure( double machNumber );

//! Compute pressure coefficient using the ACM empirical method.
/*! Computes tangent cone pressure coefficient based on the ACM empirical
 * method.
 * \param inclinationAngle Angle between wall and freestream velocity vector.
 * \param machNumber Flow Mach number.
 * \return ACM empirical surface pressure coefficient.
 */
double computeAcmEmpiricalPressureCoefficient(
        double inclinationAngle, double machNumber );

//! Compute Mach number from Prandtl-Meyer function.
/*!
 * Computes the inverse of the Prandtl-Meyer function. Currently the function
 * is limited to use with ratio of specific heats of 1.4 as it uses an
 * empirical correlation.
 * \param prandtlMeyerFunctionValue Prandyl-Meyer function value.
 * \return Mach number.
 */
double computeInversePrandtlMeyerFunction( double prandtlMeyerFunctionValue );

//! Compute ratio of post- to pre-shock pressure.
/*!
 * Computes ratio of post- to pre-shock pressure, assuming thermally and calorically perfect gas.
 * \param normalMachNumber Mach number of flow velocity normal to shock.
 * \param ratioOfSpecificHeats Ratio of specific heat at constant pressure to
 *          specific heat at constant volume.
 */
double computeShockPressureRatio( double normalMachNumber,
                                  double ratioOfSpecificHeats );

//! Compute ratio of post- to pre-shock density.
/*!
 * Compute ratio of post- to pre-shock density, assuming thermally and calorically perfect gas.
 * \param normalMachNumber Mach number of flow velocity normal to shock.
 * \param ratioOfSpecificHeats Ratio of specific heat at constant pressure to
 *         specific heat at constant volume.
 */
double computeShockDensityRatio( double normalMachNumber,
                                 double ratioOfSpecificHeats );

//! Compute ratio of post- to pre-shock temperature.
/*!
 * Computes ratio of post- to pre-shock temperature, assuming thermally and calorically perfect
 * gas.
 * \param normalMachNumber Mach number of flow velocity normal to shock.
 * \param ratioOfSpecificHeats Ratio of specific heat at constant pressure to
 *         specific heat at constant volume.
 */
double computeShockTemperatureRatio( double normalMachNumber,
                                     double ratioOfSpecificHeats );

//! Compute jump in entropy across a shock wave.
/*!
 * Compute jump in entropy across a shock wave, assuming thermally and calorically perfect gas.
 * \param normalMachNumber Mach number of flow velocity normal to shock.
 * \param ratioOfSpecificHeats ratio of specific heat at constant pressure to
 *         specific heat at constant volume.
 * \param specificGasConstant gas constant per unit mass for flow composition.
 */
double computeShockEntropyJump( double normalMachNumber, double ratioOfSpecificHeats,
                                double specificGasConstant );

//! Compute post- to pre-shock total pressure ratio.
/*!
 * Compute post- to pre-shock total pressure ratio from the entropy jump across
 * a shock wave. Assumption of thermally and calorically perfect gas is made.
 * \param normalMachNumber Mach number of flow velocity normal to shock.
 * \param ratioOfSpecificHeats ratio of specific heat at constant pressure to
 *         specific heat at constant volume.
 * \param specificGasConstant gas constant per unit mass for flow composition.
 */
double computeShockTotalPressureRatio( double normalMachNumber,
                                       double ratioOfSpecificHeats,
                                       double specificGasConstant );

//! Compute shock deflection angle.
/*!
 * Computes the flow deflection angle across a shock wave.
 * \param shockAngle Angle of shock wave w.r.t. freestream flow.
 * \param machNumber Freestream Mach number.
 * \param ratioOfSpecificHeats Ratio of specific heat at constant pressure to
 *          specific heat at constant volume.
 */
double computeShockDeflectionAngle( double shockAngle, double machNumber,
                                    double ratioOfSpecificHeats );

//! Function to compute the speed of sound in a gas
/*!
 * Function to compute the speed of sound in a gas
 * \param temperature Gas temperature
 * \param ratioOfSpecificHeats Ratio of specific heats aat constant pressure and constant volume
 * \param specificGasConstant Specific gas constant of the gas
 * \return Speed of sound in the gas.
 */
double computeSpeedOfSound( const double temperature, const double ratioOfSpecificHeats,
                            const double specificGasConstant );

//! Compute Mach number
/*!
 * Compute Mach number
 * \param speed Airspeed of object for which Mach number is to be computed.
 * \param speedOfSound Speed of sound for atmosphere position at which Mach number is to be computed.
 * \return Mach number
 */
double computeMachNumber( const double speed, const double speedOfSound );

//! Function to compute the mean free path of a particle.
/*!
 * Function to compute the mean free path of a particle from e.g. (Chapman, S. & Cowling, T. The mathematical theory of
 * nonuniform gases Cambridge University Press, 1970)
 * \param weightedAverageCollisionDiameter Weighted (using specie number density) average collision diameter of the
 * particles in the gas.
 * \param averageNumberDensity Average number density of the gas.
 * \return Mean free path of a particle in the gas.
 */
double computeMeanFreePath( const double weightedAverageCollisionDiameter, const double averageNumberDensity );

//! Function to compute the aerodynamic load experienced by a vehicle.
/*!
 * Function that computes the aerodynamic load (a.k.a. load factor) experienced by a vehicle.
 * \param airDensity Freestream air density.
 * \param airSpeed Airspeed of the vehicle.
 * \param referenceArea Reference area of the vehicle.
 * \param vehicleMass Mass of the vehicle.
 * \param aerodynamicForceCoefficients Aerodynamic force coefficients of the vehicle.
 * \return Aerodynamic load experienced by the vehicle.
 */
double computeAerodynamicLoad( const double airDensity,
                               const double airSpeed,
                               const double referenceArea,
                               const double vehicleMass,
                               const Eigen::Vector3d& aerodynamicForceCoefficients );

//! Function to compute the aerodynamic load experienced by a vehicle.
/*!
 * Function that computes the aerodynamic load (a.k.a. load factor) experienced by a vehicle.
 * \param aerodynamicAccelerationVector Aerodynamic acceleration actinv on vehicle.
 * \return Aerodynamic load experienced by the vehicle.
 */
double computeAerodynamicLoadFromAcceleration( const Eigen::Vector3d& aerodynamicAccelerationVector );

//! Funtion to compute the equilibrium heat flux experienced by a vehicle
/*!
 * Funtion to compute the equilibrium heat flux experienced by a vehicle.
 * \param heatTransferFunction Function returning the feat flux as a function of wall temperature.
 * \param wallEmmisivity Emmissivity of the wall to which heat transfer is taking place
 * \param adiabaticWallTemperature Adiabatic wall temperature (used only for initialization of root finder).
 * \return Convective heat flux experienced acording to Fay Riddell model at equilibrium wall temperature.
 */
double computeEquilibriumHeatflux( const std::function< double( const double ) > heatTransferFunction,
                                   const double wallEmmisivity,
                                   const double adiabaticWallTemperature );

//! Function to compute the heat flux experienced by a vehicle, assuming an equlibrium wall temperature.
/*!
 * Function that computes the heat flux experienced by a vehicle. This function is an implementation of the
 * Fay-Riddell formula, assuming an equlibrium wall temperature
 * \param airDensity Freestream density of the air.
 * \param airSpeed Airspeed of the vehicle.
 * \param airTemperature Freestream air temperature.
 * \param machNumber Freestream Mach number.
 * \param noseRadius Nose radius of the vehicle.
 * \param wallEmissivity Wall emissivity of the vehicle.
 * \return Convective heat flux experienced by the vehicle.
 */
double computeEquilibriumFayRiddellHeatFlux( const double airDensity,
                                             const double airSpeed,
                                             const double airTemperature,
                                             const double machNumber,
                                             const double noseRadius,
                                             const double wallEmissivity = 0.80 );

static const double FAY_RIDDEL_HEAT_FLUX_CONSTANT = 3.53E-4;

//! Function to compute the heat flux experienced by a vehicle.
/*!
 * Function that computes the heat flux experienced by a vehicle. This function is an implementation of the
 * Fay-Riddell formula.
 * \param airDensity Freestream density of the air.
 * \param airSpeed Airspeed of the vehicle.
 * \param airTemperature Freestream air temperature.
 * \param noseRadius Nose radius of the vehicle.
 * \param wallTemperature Temperature at the wall of the vehicle.
 * \return Convective heat flux experienced by the vehicle.
 */
double computeFayRiddellHeatFlux( const double airDensity,
                                  const double airSpeed,
                                  const double airTemperature,
                                  const double noseRadius,
                                  const double wallTemperature );

//! Function to compute the adiabatic wall temperature experienced by a vehicle.
/*!
 * Function that computes the adiabatic wall temperature experienced by a vehicle.
 * \param airTemperature Freestream air temperature.
 * \param machNumber Freestream Mach number.
 * \param ratioSpecificHeats Ratio of specific heats of the air.
 * \param recoveryFactor Recovery factor of flow, e.g. fraction of total enthalpy contribution from velocoty that can be
 * recovered at the wal
 * \return Adiabatic wall temperature experienced by the vehicle.
 */
double computeAdiabaticWallTemperature(
        const double airTemperature, const double machNumber, const double ratioSpecificHeats = 1.4,
        const double recoveryFactor = 0.845 );

} // namespace aerodynamics

} // namespace tudat

#endif // TUDAT_AERODYNAMICS_H

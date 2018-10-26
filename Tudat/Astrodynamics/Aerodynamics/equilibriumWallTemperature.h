/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_EQUILIBRIUMWALLTEMPERATURE_H
#define TUDAT_EQUILIBRIUMWALLTEMPERATURE_H

#include <functional>

#include "Tudat/Mathematics/BasicMathematics/basicFunction.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/basicElectroMagnetism.h"

namespace tudat
{

namespace aerodynamics
{

//! Function that is used to compute the net heat flux from a given heat input and wall temperature
class EquilibriumTemperatureFunction: public basic_mathematics::Function<double,double>
{
public:
    //! Constructor.
    /*!
     * Constructor
     * \param heatTransferFunction Function returning the feat flux as a function of wall temperature.
     * \param wallEmissivity Emmissivity of the wall to which heat transfer is taking place
     * \param adiabaticWallTemperature Adiabatic wall temperature
     */
    EquilibriumTemperatureFunction(
            const std::function< double( const double ) > heatTransferFunction,
            const double wallEmissivity,
            double adiabaticWallTemperature ):
       heatTransferFunction_( heatTransferFunction ), wallEmissivity_( wallEmissivity ),
       adiabaticWallTemperature_( adiabaticWallTemperature ){ }

    //! Destructor.
    ~EquilibriumTemperatureFunction(){}

    //! Compute net heat flux at given wall temperature
    /*!
     * Compute net heat flux at given wall temperature
     * \param currentWallTemperature Wall temperature to be used for computation of input and output of heat.
     * \return Net heat input to wall
     */
    double evaluate( const double currentWallTemperature )
    {
        return heatTransferFunction_( currentWallTemperature )
                - wallEmissivity_ * electro_magnetism::computeBlackbodyRadiationIntensity( currentWallTemperature );
    }

    //! Compute first derivative of net heat flux at given wall temperature (FUNCTION NOT IMPLEMENTED)
    double computeDerivative( const unsigned int order, const double independentVariable )
    {
        throw std::runtime_error( "Error, derivative of heat flux not defined" );
        return TUDAT_NAN;
    }

    //! Compute first derivative of net heat flux at given wall temperature (FUNCTION NOT IMPLEMENTED)
    double computeDefiniteIntegral( const unsigned int order, const double lowerBound, const double upperbound )
    {
        throw std::runtime_error( "Error, integrall of heat flux not defined" );
        return TUDAT_NAN;
    }

    //! Function to retrieve the lower bound of the wall temperature.
    double getLowerBound( ) { return 0.0; }

    //! Function to retrieve the upper bound of the wall temperature.
    double getUpperBound( ) { return adiabaticWallTemperature_; }

    //! Function to retrieve the initial guess of the wall temperature.
    double getInitialGuess( ) { return adiabaticWallTemperature_ * 0.01; }

protected:

private:
    //! Function that returns the heat input as a function of wall temperature.
    std::function< double( const double ) > heatTransferFunction_;

    //! Constant wall emissivity.
    const double wallEmissivity_;

    //! Constant adiabatic wall temperature.
    double adiabaticWallTemperature_;
};

//! Function to compute the equilibrium wall temperature from the heat input and emmisivity
/*!
 * Function to compute the equilibrium wall temperature from the heat input and emmisivity. This function calls a root-
 * finder to determine the wall temperature where the heat input equals the radiative output.
 * \param heatTransferFunction Function that returns the heat input as a function of current temperature
 * \param wallEmmisivity Value of the material emmisivity
 * \param adiabaticWallTemperature Adiabatic wall temperature. This variables is only used as an upper bound, and to
 * generate an initial giess for the equilibrium temperature
 * \return Wall temperature at which input and output of heat are in equilibrium.
 */
double computeEquilibiumWallTemperature(
        const std::function< double( const double ) > heatTransferFunction,
        const double wallEmmisivity,
        const double adiabaticWallTemperature );

} //namespace_aerodynamics

} //namespace_tudat

#endif //TUDAT_EQUILIBRIUMWALLTEMPERATURE_H

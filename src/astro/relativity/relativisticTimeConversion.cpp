/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#include <iostream>
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/relativity/relativisticTimeConversion.h"
#include <iomanip>

namespace tudat
{

namespace relativity
{

//! Function to compute proper-time rate w.r.t. coordinate time, minus 1.0, from a speed and scalar potential
double calculateFirstCentralBodyProperTimeRateDifference(
        const double computationPointSpeed,
        double gravitationalScalarPotential,
        const double equivalencePrincipleLpiViolationParameter )
{
    return ( -( 0.5 * computationPointSpeed * computationPointSpeed + ( 1.0 + equivalencePrincipleLpiViolationParameter ) *
                gravitationalScalarPotential )  ) * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
}

//! Function to compute proper-time rate w.r.t. coordinate time, minus 1.0, from a speed and scalar potential
double calculateFirstCentralBodyProperTimeRateDifference(
    const Eigen::Vector6d& computationPointState,
    const std::vector< Eigen::Vector6d >& perturbedInertialStates,
    const std::vector< double >& centralBodyGravitationalParameters,
    const double equivalencePrincipleLpiViolationParameter )
{
    double gravitationalScalarPotential = 0.0;
    for( unsigned int i = 0; i < perturbedInertialStates.size( ); i++ )
    {
        gravitationalScalarPotential +=
            centralBodyGravitationalParameters.at( i ) /
                ( perturbedInertialStates.at( i ).segment( 0, 3 ) - computationPointState.segment( 0, 3 ) ).norm( );
    }

    return calculateFirstCentralBodyProperTimeRateDifference(
        computationPointState.segment( 3, 3 ).norm( ), gravitationalScalarPotential, equivalencePrincipleLpiViolationParameter );
}

}

}


/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/aerodynamics/equilibriumWallTemperature.h"
#include "tudat/math/root_finders/rootFinder.h"
#include "tudat/math/root_finders/terminationConditions.h"
#include "tudat/math/root_finders/bisection.h"
#include "tudat/math/root_finders/secantRootFinder.h"

namespace tudat
{

namespace aerodynamics
{
//! Function to compute the equilibrium wall temperature from the heat input and emmisivity
double computeEquilibiumWallTemperature( const std::function< double( const double ) > heatTransferFunction,
                                         const double wallEmmisivity,
                                         const double adiabaticWallTemperature )
{
    // Create the object that contains the function who's root needs to be found.
    std::shared_ptr< EquilibriumTemperatureFunction > equilibriumTemperatureFunction
            = std::make_shared< EquilibriumTemperatureFunction >(
                heatTransferFunction, wallEmmisivity, adiabaticWallTemperature  );

    // Compute wall temperature, first try secant method, use bisection as backup if secany unsuccesfull.
    double wallTemperature = TUDAT_NAN;
    try
    {
        root_finders::SecantRootFinder< >::TerminationFunction terminationConditionFunction =
                std::bind( &root_finders::RootRelativeToleranceTerminationCondition< double >::
                             checkTerminationCondition,
                             std::make_shared< root_finders::
                             RootRelativeToleranceTerminationCondition< double > >(
                                 ), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5 );
        root_finders::SecantRootFinder< > secant( terminationConditionFunction );
        wallTemperature = secant.execute(
                    equilibriumTemperatureFunction, equilibriumTemperatureFunction->getInitialGuess( ) );
    }
    catch( std::runtime_error const& )

    {
        try
        {
        root_finders::Bisection< >::TerminationFunction terminationConditionFunction =
                std::bind( &root_finders::RootRelativeToleranceTerminationCondition< double >::
                             checkTerminationCondition,
                             std::make_shared< root_finders::
                             RootRelativeToleranceTerminationCondition< double > >(
                                 ), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5 );
        root_finders::Bisection< > bisection( terminationConditionFunction );
        wallTemperature = bisection.execute(
                    equilibriumTemperatureFunction, equilibriumTemperatureFunction->getInitialGuess( ) );
        }
        catch( std::runtime_error const& )

        {
            throw std::runtime_error( "Error, could not find equilibrium wall temperature" );
        }

    }

    return wallTemperature;
}

} //namespace_aerodynamics

} //namespace_tudat


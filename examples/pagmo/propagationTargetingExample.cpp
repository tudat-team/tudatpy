/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


/* A satellite is going to be launched on an elliptical orbit with min and max altitudes
respectively 180 km and 40000 km with inclination i=0 (and argument of perigee omega = 0).
A fixed target is set on the same plane at an altitude of 35000 km, fixed latitude of 30 deg.
What is the value of the RAAN for which the satellite achieves a minimum
approach distance from the target?*/

#include <pagmo/problem.hpp>
#include <pagmo/algorithms/sade.hpp>
#include <pagmo/algorithms/de1220.hpp>
#include <pagmo/algorithms/de.hpp>
#include <pagmo/algorithms/simulated_annealing.hpp>
#include <pagmo/io.hpp>
#include <pagmo/archipelago.hpp>

#include "Problems/propagationTargeting.h"
#include "Problems/applicationOutput.h"
#include "Problems/saveOptimizationResults.h"

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/io/basicInputOutput.h"

using namespace pagmo;
using namespace tudat_pagmo_applications;
using namespace tudat;

int main( )
{
    bool performGridSearch = false;

    //Set seed for reproducible results
    pagmo::random_device::set_seed(255);

    //Load spice kernels
    tudat::spice_interface::loadStandardSpiceKernels( );

    // Create object to compute the problem fitness; no perturbations
    double altitudeOfPerigee = 180000.0;
    double altitudeOfApogee = 40000000.0;
    double altitudeOfTarget = 35000000.0;
    double longitudeOfTarget = 30.0; // In degrees

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                                          propagators::altitude_dependent_variable, "Satellite", "Earth" ) );

    // Create object with list of dependent variables.
    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< propagators::DependentVariableSaveSettings >( dependentVariablesList );


    PropagationTargetingProblem targetingProblem( altitudeOfPerigee, altitudeOfApogee, altitudeOfTarget,
                                                  longitudeOfTarget, dependentVariablesToSave, false );

    problem prob{ targetingProblem };

    // Perform Grid Search and write results to file
    if( performGridSearch )
    {
        createGridSearch( prob, {{0.0, 0.0}, {360.0, 180.0}}, { 100, 50 }, "propagationTargetingGridSearch_" );
    }
    // Instantiate a pagmo algorithm
    algorithm algo{de1220( )};

    // Create an island with 128 individuals
    pagmo::population::size_type populationSize = 128;
    island isl{algo, prob, populationSize};

    // Evolve for 25 generations
    for( int i = 0; i < 25; i++ )
    {
        isl.evolve( );
        while( isl.status( ) != pagmo::evolve_status::idle &&
               isl.status( ) != pagmo::evolve_status::idle_error )
        {
            isl.wait( );
        }
        isl.wait_check( ); // Raises errors

        // Write current iteration results to file
        printPopulationToFile( isl.get_population( ).get_x( ), "targetingPropagation_" + std::to_string( i ) + "_" + std::to_string( i ) , false );
        printPopulationToFile( isl.get_population( ).get_f( ), "targetingPropagation_" + std::to_string( i ) + "_" + std::to_string( i ) , true );

        std::cout<<i<<std::endl;
    }

    // Retrieve final Cartesian states for population in last generation, and save final states to a file.
    std::vector<std::vector< double > > decisionVariables = isl.get_population( ).get_x( );
    std::map< int, Eigen::VectorXd > finalStates;
    for( unsigned int i = 0; i < decisionVariables.size( ); i++ )
    {
        targetingProblem.fitness( decisionVariables.at( i ) );
        finalStates[ i ] = targetingProblem.getPreviousFinalState( );
    }
    tudat::input_output::writeDataMapToTextFile(
                finalStates, "targetingFinalStates.dat", tudat_pagmo_applications::getOutputPath( ) );

    // Retrieve final values of dependent variables for population in last generation, and save final dependent variables values to a file.
    std::map< int, Eigen::VectorXd > dependentVariablesFinalValues;
    for( unsigned int i = 0; i < decisionVariables.size( ); i++ )
    {
        targetingProblem.fitness( decisionVariables.at( i ) );
        dependentVariablesFinalValues[ i ] = targetingProblem.getPreviousDependentVariablesFinalValues();
    }
    tudat::input_output::writeDataMapToTextFile(
                dependentVariablesFinalValues, "targetingDependentVariablesFinalVariables.dat", tudat_pagmo_applications::getOutputPath( ) );


    // Create object to compute the problem fitness; with perturbations
    problem prob_pert{PropagationTargetingProblem( altitudeOfPerigee, altitudeOfApogee, altitudeOfTarget,
                                                   longitudeOfTarget, dependentVariablesToSave, true ) };


    // Instantiate a pagmo algorithm for the new problem.
    algorithm algo_pert{de1220( )};

    // Create an empty population for perturbed problem
    population population_pert = population( prob_pert, 0 );

    // Retrieve population of unperturbed problem, and instantiate population of perturbed problem
    std::vector<vector_double> original_population = isl.get_population( ).get_x( );
    for( unsigned int k = 0; k < populationSize; k++ )
    {
        population_pert.push_back( original_population.at( k ) );
    }

    // Create island for perturbed problem
    island isl_pert{algo_pert, population_pert};

    // Perform Grid Search for perturbed priblem and write results to file
    if( performGridSearch )
    {
        createGridSearch( isl_pert.get_population( ).get_problem( ), {{0.0, 0.0}, {360.0, 180.0}}, { 100, 50 }, "propagationTargetingGridSearch_pert" );
    }

    // Write original (unevolved) population to file
    printPopulationToFile( isl_pert.get_population( ).get_x( ), "targetingPropagation_pert_orig" , false );
    printPopulationToFile( isl_pert.get_population( ).get_f( ), "targetingPropagation_pert_orig" , true );


    // Evolve for 4 generations
    for( int i = 0; i < 4; i++ )
    {
        isl_pert.evolve( );
        while( isl_pert.status( ) != pagmo::evolve_status::idle &&
               isl_pert.status( ) != pagmo::evolve_status::idle_error )
        {
            isl_pert.wait( );
        }
        isl_pert.wait_check( ); // Raises errors

        // Write current iteration results to file
        printPopulationToFile( isl_pert.get_population( ).get_x( ), "targetingPropagation_pert_" + std::to_string( i ) + "_" + std::to_string( i ) , false );
        printPopulationToFile( isl_pert.get_population( ).get_f( ), "targetingPropagation_pert_" + std::to_string( i ) + "_" + std::to_string( i ) , true );

        std::cout<<i<<std::endl;
    }
}


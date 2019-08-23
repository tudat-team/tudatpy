/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */


#include "Tudat/Astrodynamics/LowThrustDirectMethods/hybridMethod.h"
#include "Tudat/Astrodynamics/LowThrustDirectMethods/hybridOptimisationSetup.h"
#include "Tudat/Astrodynamics/LowThrustDirectMethods/hybridMethodLeg.h"

#include "tudatExampleApplications/libraryExamples/PaGMOEx/Problems/applicationOutput.h"
#include "tudatExampleApplications/libraryExamples/PaGMOEx/Problems/getAlgorithm.h"
#include "tudatExampleApplications/libraryExamples/PaGMOEx/Problems/saveOptimizationResults.h"
#include "pagmo/problems/unconstrain.hpp"
#include "pagmo/algorithms/compass_search.hpp"


namespace tudat
{
namespace low_thrust_direct_methods
{


//! Perform optimisation.
std::pair< std::vector< double >, std::vector< double > > HybridMethod::performOptimisation( )
{
    using namespace tudat_pagmo_applications;

    //Set seed for reproducible results
    pagmo::random_device::set_seed( 456 );

    // Create object to compute the problem fitness
    problem prob{ HybridMethodProblem( stateAtDeparture_, stateAtArrival_, maximumThrust_, specificImpulseFunction_, numberOfRevolutions_,
                                       timeOfFlight_, bodyMap_, bodyToPropagate_, centralBody_, integratorSettings_, relativeToleranceConstraints_,
                                       optimiseTimeOfFlight_, timeOfFlightBounds_ )};

    std::vector< double > constraintsTolerance;
    for ( unsigned int i = 0 ; i < ( prob.get_nec() + prob.get_nic() ) ; i++ )
    {
        constraintsTolerance.push_back( 1.0e-3 );
    }
    prob.set_c_tol( constraintsTolerance );


//    unconstrain unconstrainedProb{ prob, "ignore_o" };
//    population pop{ unconstrainedProb, 10 };

    algorithm algo = optimisationAlgorithm_;
//    algoUnconstrained.set_verbosity( 10 );
//    algoUnconstrained.evolve( pop );

    island island{ algo, prob, 50 /*5000*/ };

    // Evolve for 10 generations
    for( int i = 0 ; i < 2 /*600*/ ; i++ )
    {
        island.evolve( );
        while( island.status( ) != pagmo::evolve_status::idle &&
               island.status( ) != pagmo::evolve_status::idle_error )
        {
            island.wait( );
        }
        island.wait_check( ); // Raises errors

        // Write current iteration results to file
        printPopulationToFile( island.get_population( ).get_x( ), "testHybridMethod_generation_" + std::to_string( i ) , false );
        printPopulationToFile( island.get_population( ).get_f( ), "testHybridMethod_generation_" + std::to_string( i ) , true );

        std::vector< double > championFitness = island.get_population().champion_f();
        for ( int i = 0 ; i < championFitness.size() ; i++ )
        {
            std::cout << "champion fitness: " << championFitness[ i ] << "\n\n";
        }
        std::cout << "TEST" << "\n\n";
        std::cout<< "current generation: " << i << std::endl;
    }
    std::cout << "TEST" << "\n\n";

    std::vector< double > championFitness = island.get_population().champion_f();
    std::vector< double > championDesignVariables = island.get_population().champion_x();
    for ( int i = 0 ; i < championFitness.size() ; i++ )
    {
        std::cout << "champion fitness: " << championFitness[ i ] << "\n\n";
    }
    for ( int i = 0 ; i < championDesignVariables.size() ; i++ )
    {
        std::cout << "champion design variables: " << championDesignVariables[ i ] << "\n\n";
    }

    std::vector< double > championFitnessConstrainedPb = prob.fitness( championDesignVariables );
    for ( int i = 0 ; i < championFitnessConstrainedPb.size() ; i++ )
    {
        std::cout << "champion fitness constrained problem: " << championFitnessConstrainedPb[ i ] << "\n\n";
    }

    std::cout << "champion fitness: " << championFitness[ 0 ] << "\n\n";



//        // Create an island with 1024 individuals
//        island isl{ algo, prob, 100}; //1024};

//        // Evolve for 100 generations
//        for( int i = 0 ; i < 1 ; i++ ) //300 ; i++) //100; i++ )
//        {
//            isl.evolve( );
//            while( isl.status( ) != pagmo::evolve_status::idle &&
//                   isl.status( ) != pagmo::evolve_status::idle_error )
//            {
//                isl.wait( );
//            }
//            isl.wait_check( ); // Raises errors

//            // Write current iteration results to file
//            printPopulationToFile( isl.get_population( ).get_x( ), "testSimsFlanagan_generation_" + std::to_string( i ) , false );
//            printPopulationToFile( isl.get_population( ).get_f( ), "testSimsFlanagan_generation_" + std::to_string( i ) , true );

//            std::cout<< "current generation: " << i << std::endl;
//        }

    std::pair< std::vector< double >, std::vector< double > > output;
    output.first = championFitness;
    output.second = championDesignVariables;

    championFitness_ = championFitness;
    championDesignVariables_ = championDesignVariables;

    return output;

}


} // namespace low_thrust_direct_methods
} // namespace tudat

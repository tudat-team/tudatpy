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
#include <fstream>

#include <boost/filesystem.hpp>

#include "pagmo/algorithms/de1220.hpp"
#include "pagmo/algorithms/sade.hpp"
#include "pagmo/algorithms/cmaes.hpp"
#include "pagmo/algorithms/moead.hpp"
#include "pagmo/algorithms/nsga2.hpp"
#include "pagmo/algorithms/ihs.hpp"
#include "Problems/earthMarsTransfer.h"
#include "Problems/applicationOutput.h"
#include "Problems/getAlgorithm.h"
#include "Problems/saveOptimizationResults.h"

using namespace tudat_pagmo_applications;

int main( )
{

    //Set seed for reproducible results
    pagmo::random_device::set_seed( 123 );

    // We have two decision variables each with a lower and upper bound, create a vector of vectors that will contain these.
    std::vector< std::vector< double > > bounds( 2, std::vector< double >( 2, 0.0 ) );

    // Define bounds: Search between 2020 and 2025 for flight duration between 200 and 1000 days.
    bounds[ 0 ][ 0 ] = 2458849.5;
    bounds[ 1 ][ 0 ] = 2460676.5;
    bounds[ 0 ][ 1 ] = 200;
    bounds[ 1 ][ 1 ] = 1000;

    // Create object to compute the problem fitness
    problem prob{EarthMarsTransfer( bounds, true )};

    // Solve problem using 3 different optimizers
    for( unsigned int j = 0; j < 3; j++ )
    {
        // Retrieve MO algorithm
        algorithm algo{getMultiObjectiveAlgorithm( j )};

        // Create an island with 1024 individuals
        island isl{algo, prob, 1024};

        // Evolve for 25 generations
        for( int i = 0 ; i < 100; i++ )
        {
            isl.evolve( );
            while( isl.status( ) != pagmo::evolve_status::idle &&
                   isl.status( ) != pagmo::evolve_status::idle_error )
            {
                isl.wait( );
            }
            isl.wait_check( ); // Raises errors

            // Write current iteration results to file
            printPopulationToFile( isl.get_population( ).get_x( ), "mo_EarthMars_" + std::to_string( i ) + "_" + std::to_string( j ), false );
            printPopulationToFile( isl.get_population( ).get_f( ), "mo_EarthMars_" + std::to_string( i ) + "_" + std::to_string( j ), true );

            std::cout<<i<<" "<<j<<std::endl;


        }
    }

    return 0;

}

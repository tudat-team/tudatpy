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
#include "pagmo/algorithms/simulated_annealing.hpp"
#include "pagmo/algorithms/cmaes.hpp"
#include "pagmo/algorithms/de.hpp"
#include "pagmo/algorithms/nsga2.hpp"

#include "Problems/multipleGravityAssist.h"
#include "Problems/applicationOutput.h"
#include "Problems/getAlgorithm.h"
#include "Problems/saveOptimizationResults.h"

//! Execute  main
int main( )
{
    //Set seed for reproducible results
    pagmo::random_device::set_seed( 123456789 );

    // We have five decision variables each with a lower and upper
    // bound, create a vector of vectors that will contain these.
    int numberOfParameters = 5;
    std::vector< std::vector< double > > bounds( 2, std::vector< double >( numberOfParameters, 0.0 ) );

    // Define search bounds: first parameter is start date, following parameters are leg durations
    bounds[ 0 ][ 0 ] = 7304.5; //MJD2000
    bounds[ 1 ][ 0 ] = 7304.5 + 10 * 365; //MJD2000
    bounds[ 0 ][ 1 ] = 200;
    bounds[ 1 ][ 1 ] = 500;
    bounds[ 0 ][ 2 ] = 50;
    bounds[ 1 ][ 2 ] = 300;
    bounds[ 0 ][ 3 ] = 50;
    bounds[ 1 ][ 3 ] = 300;
    bounds[ 0 ][ 4 ] = 500;
    bounds[ 1 ][ 4 ] = 3000;

    // Define the problem: EVEEJ flyby sequence
    std::vector< int > flybySequence;
    flybySequence.push_back( 3 );
    flybySequence.push_back( 2 );
    flybySequence.push_back( 3 );
    flybySequence.push_back( 3 );
    flybySequence.push_back( 5 );

    // Create object to compute the problem fitness
    problem prob{ MultipleGravityAssist( bounds, flybySequence, true ) };

    // Select NSGA2 algorithm for priblem
    algorithm algo{nsga2( )};

    // Create an island with 1000 individuals
    island isl{algo, prob, 1000 };

    // Evolve for 512 generations
    for( int i = 0 ; i < 512; i++ )
    {

        isl.evolve( );
        while( isl.status( ) != pagmo::evolve_status::idle &&
               isl.status( ) != pagmo::evolve_status::idle_error )
        {
            isl.wait( );
        }

        isl.wait_check( ); // Raises errors

        // Write current iteration results to file
        printPopulationToFile( isl.get_population( ).get_x( ), "mo_mga_EVEEJ_" + std::to_string( i ), false );
        printPopulationToFile( isl.get_population( ).get_f( ), "mo_mga_EVEEJ_" + std::to_string( i ), true );
        std::cout<<i<<std::endl;
    }

    return 0;

}

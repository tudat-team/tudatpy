#include <iostream>
#include "pagmo/algorithms/de1220.hpp"
#include "pagmo/algorithms/pso.hpp"
#include "pagmo/algorithms/de.hpp"
#include "pagmo/algorithms/sga.hpp"
#include "pagmo/algorithms/sade.hpp"
#include "pagmo/island.hpp"
#include "pagmo/problem.hpp"
#include "Problems/himmelblau.h"
#include "Problems/applicationOutput.h"
#include "Problems/saveOptimizationResults.h"

int main( )
{
    using namespace tudat_pagmo_applications;

    //Set seed for reproducible results
    pagmo::random_device::set_seed( 12345 );

    // Define Himmelblau function (range [-5,5]x[-5,5])
    pagmo::problem prob{ HimmelblauFunction( -5, 5, -5, 5) };

    // Perform grid saerch
    //createGridSearch( prob, { {- 5.0, -5.0 },{ 5.0, 5.0 } }, { 1000, 1000 }, "himmelBlauGridSearch" );

    // Solve using DE algorithm
    pagmo::algorithm algo{ pagmo::de( ) };

    // Create island with 1000 individuals
    pagmo::island isl = pagmo::island{ algo, prob, 1000 };

    // Evolve for 1000 generations
    for( int i = 1; i <= 100; i++ )
    {
        isl.evolve( );
        while( isl.status( ) != pagmo::evolve_status::idle &&
               isl.status( ) != pagmo::evolve_status::idle_error )
        {
            isl.wait( );
        }
        isl.wait_check( ); // Raises errors

        printPopulationToFile( isl.get_population( ).get_x( ), "himmelblau_" + std::to_string( i ) , false );
        printPopulationToFile( isl.get_population( ).get_f( ), "himmelblau_" +  std::to_string( i ) , true );

        // Print current optimum to console
        std::cout << "Minimum: " <<i<<" "<<std::setprecision( 16 ) <<"f= "<< isl.get_population().champion_f()[0] <<", x="<<
                     isl.get_population().champion_x()[0] <<" y="<<isl.get_population().champion_x()[1] <<std::endl;

    }


    return 0;

}

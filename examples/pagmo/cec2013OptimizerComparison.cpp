#include <iostream>

#include <boost/filesystem.hpp>

#include "pagmo/problems/cec2013.hpp"
#include "pagmo/island.hpp"
#include "pagmo/problem.hpp"

#include "Problems/himmelblau.h"
#include "Problems/getAlgorithm.h"
#include "Problems/saveOptimizationResults.h"
#include "Problems/applicationOutput.h"

int main( )
{
    // Set random seed
    pagmo::random_device::set_seed( 12345 );

    // Define run settings
    unsigned int numberOfDimensionCases = 3;
    unsigned int numberOfOptimizers = 11;
    unsigned int numberOfProblems = 28;
    unsigned int numberOfPopulationCases = 3;

    // Create double vector of Matrices that contain optimal solutions
    std::vector< std::vector< Eigen::MatrixXd > > optima;
    optima.resize( numberOfDimensionCases );
    for( unsigned int i = 0; i < numberOfDimensionCases; i++ )
    {
        optima[ i ].resize( numberOfPopulationCases );
        for( unsigned int j = 0; j < numberOfPopulationCases; j++ )
        {
            optima[ i ][ j ] = Eigen::MatrixXd( numberOfProblems, numberOfOptimizers );
        }
    }

    // Perform optimization for each optimizer
    for( unsigned int i = 0; i < numberOfOptimizers; i++ )
    {
        // Perform optimization for each test problem
        for( unsigned int j = 0; j < numberOfProblems; j++ )
        {
            // Perform optimization for each dimensionality
            for( unsigned int k = 0; k < numberOfDimensionCases; k++ )
            {
                int numberOfDimensions = 0;
                if( k == 0 )
                {
                    numberOfDimensions = 2;
                }
                else if( k == 1 )
                {
                    numberOfDimensions = 5;
                }
                else if( k == 2 )
                {
                    numberOfDimensions = 10;
                }

                // Perform optimization for each population size/number of generations
                for( unsigned int l = 0; l < numberOfPopulationCases; l++ )
                {
                    pagmo::population::size_type populationSize = 0;
                    int numberOfGenerations = 0;
                    if( l == 0 )
                    {
                        populationSize = 16;
                        numberOfGenerations = 64;
                    }
                    else if( l == 1 )
                    {
                        populationSize = 32;
                        numberOfGenerations = 32;
                    }
                    else if( l == 2 )
                    {
                        populationSize = 64;
                        numberOfGenerations = 16;
                    }

                    // Create test problem
                    pagmo::problem prob{ pagmo::cec2013( j + 1, numberOfDimensions ) };

                    // Create optimization algorithm
                    pagmo::algorithm algo = getAlgorithm( i );

                    // Create island with given population size
                    pagmo::island isl = pagmo::island{ algo, prob, populationSize };

                    // Evolve for required number of generations
                    for( int j = 1; j <= numberOfGenerations; j++ )
                    {
                        isl.evolve( );
                        while( isl.status( ) != pagmo::evolve_status::idle &&
                               isl.status( ) != pagmo::evolve_status::idle_error )
                        {
                            isl.wait( );
                        }
                        isl.wait_check( ); // Raises errors
                    }

                    // Save optimal results for current optimization.
                    optima[ k ][ l ]( j, i ) = isl.get_population().champion_f()[0];
                    std::cout << "Minimum: " << i<<" "<<" "<<j<<" "<<isl.get_population().champion_f()[0] << std::endl;
                }
            }
            std::cout<<std::endl;
        }
    }

    // Write all results to files.
    for( unsigned int i = 0; i < numberOfDimensionCases; i++ )
    {
        for( unsigned int j = 0; j < numberOfPopulationCases; j++ )
        {
            tudat::input_output::writeMatrixToFile( optima[ i ][ j ], "cec2013Optima_" + std::to_string( i ) + "_" +
                                                    std::to_string( j ) + ".dat", 16,
                                                    tudat_pagmo_applications::getOutputPath( ) );
        }
    }

    return 0;
}

#include <iostream>

#include <boost/filesystem.hpp>

#include "applicationOutput.h"

#include "tudat/io/basicInputOutput.h"

void printPopulationToFile( const std::vector< std::vector< double > >& population,
                            const std::string fileSuffix,
                            const bool isFitness )
{

    Eigen::MatrixXd matrixToPrint( population.size( ), population.at( 0 ).size( ) );
    for( unsigned int i = 0; i < population.size( ); i++ )
    {
        for( unsigned int j = 0; j < population.at( 0 ).size( ); j++ )
        {
            matrixToPrint( i, j ) = population.at( i ).at( j );
        }
    }

    if( !isFitness )
    {
        tudat::input_output::writeMatrixToFile( matrixToPrint, "population_" + fileSuffix + ".dat", 16,
                                                tudat_pagmo_applications::getOutputPath( ) );
    }
    else
    {
        tudat::input_output::writeMatrixToFile( matrixToPrint, "fitness_" + fileSuffix + ".dat", 16,
                                                tudat_pagmo_applications::getOutputPath( )  );
    }
}

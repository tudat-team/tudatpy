#include <iostream>

#include <boost/filesystem.hpp>

#include "pagmo/algorithms/de1220.hpp"
#include "pagmo/algorithms/pso.hpp"
#include "pagmo/algorithms/de.hpp"
#include "pagmo/algorithms/bee_colony.hpp"
#include "pagmo/algorithms/sga.hpp"
#include "pagmo/algorithms/sea.hpp"
#include "pagmo/algorithms/sade.hpp"
#include "pagmo/algorithms/simulated_annealing.hpp"
#include "pagmo/algorithms/bee_colony.hpp"
#include "pagmo/algorithms/cmaes.hpp"
#include "pagmo/algorithms/nlopt.hpp"
#include "pagmo/algorithms/nsga2.hpp"
#include "pagmo/algorithms/moead.hpp"
#include "pagmo/problems/griewank.hpp"
#include "pagmo/problems/schwefel.hpp"
#include "pagmo/problems/rosenbrock.hpp"
#include "pagmo/problems/cec2013.hpp"
#include "pagmo/problems/zdt.hpp"
#include "pagmo/island.hpp"
#include "pagmo/problem.hpp"
#include "Problems/himmelblau.h"

template< typename OutputStream, typename ScalarType,
          int NumberOfRows, int NumberOfColumns, int Options, int MaximumRows, int MaximumCols >
void writeValueToStream( OutputStream& stream, const Eigen::Matrix< ScalarType,
                         NumberOfRows, NumberOfColumns, Options,
                         MaximumRows, MaximumCols >& value,
                         const int precision, const std::string& delimiter,
                         const bool endLineAfterRow  = 0 )
{
    for ( int i = 0; i < value.rows( ); i++ )
    {
        for ( int j = 0; j < value.cols( ); j++ )
        {
            stream << delimiter << " "
                   << std::setprecision( precision ) << std::left
                   << std::setw( precision + 1 )
                   << value( i, j );


        }
        if( endLineAfterRow )
        {
            stream << std::endl;
        }
    }
    stream << std::endl;
}

template< typename ScalarType, int NumberOfRows, int NumberOfColumns >
void writeMatrixToFile( Eigen::Matrix< ScalarType, NumberOfRows, NumberOfColumns > matrixToWrite,
                        const std::string& outputFilename,
                        const int precisionOfMatrixEntries = 16,
                        const boost::filesystem::path& outputDirectory =
        "/home/dominic/Software/optimizationBundle/tudatBundle/tudatExampleApplications/libraryExamples/Pagmo2/bin/applications/",
                        const std::string& delimiter = "\t",
                        const std::string& header = "" )
{
    // Check if output directory exists; create it if it doesn't.
    if ( !boost::filesystem::exists( outputDirectory ) )
    {
        boost::filesystem::create_directories( outputDirectory );
    }

    // Open output file.
    std::string outputDirectoryAndFilename = outputDirectory.string( ) + "/" + outputFilename;
    std::ofstream outputFile_( outputDirectoryAndFilename.c_str( ) );

    // Write header
    outputFile_ << header;

    writeValueToStream( outputFile_, matrixToWrite, precisionOfMatrixEntries,
                        delimiter, true );

    outputFile_.close( );
}

pagmo::algorithm getAlgorithm( const int index )
{
    switch( index )
    {
    case 0:
    {
        pagmo::algorithm algo{ pagmo::nsga2( ) };
        return algo;
        break;
    }
    case 1:
    {
        pagmo::algorithm algo{ pagmo::moead( ) };
        return algo;
        break;
    }
    case 2:
    {
        pagmo::algorithm algo{ pagmo::ihs( ) };
        return algo;
        break;
    }
    }
}

void printPopulationToFile( const int problemIndex, const int iterationIndex,
                            const std::vector< pagmo::vector_double >& population,
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
        writeMatrixToFile( matrixToPrint, "population_mo_" + std::to_string( problemIndex ) + "_" + std::to_string( iterationIndex ) + ".dat" );
    }
    else
    {
        writeMatrixToFile( matrixToPrint, "fitness_mo_" + std::to_string( problemIndex ) + "_" + std::to_string( iterationIndex ) + ".dat" );
    }
}

int main( )
{
    pagmo::random_device::set_seed( 12345 );

    for( unsigned int i = 0; i < 3; i++ )
    {
        pagmo::problem prob{ pagmo::zdt( 3, 2 ) };//my_problem( 0, 5, 0, 5) };

        pagmo::algorithm algo = getAlgorithm( i );

        pagmo::island isl = pagmo::island{ algo, prob, 64 };

        for( int j = 1; j <= 64; j++ )
        {

            isl.evolve( );

            //isl.get_population( ).get_f()

            printPopulationToFile( i, j, isl.get_population( ).get_x( ), false );
            printPopulationToFile( i, j, isl.get_population( ).get_f( ), true );

            //std::cout << "Best x: " << isl.get_population().champion_x()[0] << std::endl;
            //std::cout << "Best y: " << isl.get_population().champion_x()[1] << std::endl;
            while( isl.status()!=pagmo::evolve_status::idle )
                isl.wait();

        }
        std::cout<<i<<" done"<<std::endl;
    }
    return 0;
}

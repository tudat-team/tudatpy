#ifndef TUDAT_PAGMO_APPLICATIONOUTPUT_H
#define TUDAT_PAGMO_APPLICATIONOUTPUT_H

#include <pagmo/problem.hpp>

#include "tudat/io/basicInputOutput.h"
#include "tudat/basics/utilities.h"

namespace tudat_pagmo_applications
{

//! Get path for output directory.
static inline std::string getOutputPath(
        const std::string& extraDirectory = "" )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string reducedPath = filePath_.substr( 0, filePath_.length( ) -
                                std::string( "Problems/applicationOutput.h" ).length( ) );
    std::string outputPath = reducedPath + "SimulationOutput/";
    if( extraDirectory != "" )
    {
        outputPath += extraDirectory;
    }

    if( outputPath.at( outputPath.size( ) - 1 ) != '/' )
    {
        outputPath += "/";
    }

    return outputPath;
}

void createGridSearch(
        pagmo::problem& problem,
        const std::vector< std::vector< double > >& bounds,
        const std::vector< int > numberOfPoints,
        const std::string& fileName )
{
    if( bounds.at( 0 ).size( ) != 2 )
    {
        std::cerr<<"Warning when plotting grid search, size of problem does not equal 2"<<std::endl;
    }
    Eigen::MatrixXd gridSearch = Eigen::MatrixXd( numberOfPoints.at( 0 ), numberOfPoints.at( 1 ) );

    double xSpacing = ( bounds[ 1 ][ 0 ] - bounds[ 0 ][ 0 ] ) / static_cast< double >( numberOfPoints.at( 0 ) - 1 );
    double ySpacing = ( bounds[ 1 ][ 1 ] - bounds[ 0 ][ 1 ] ) / static_cast< double >( numberOfPoints.at( 1 ) - 1 );

    std::vector< double > decisionVector;
    decisionVector.resize( 2 );

    std::vector< double > xDataPoints;
    for( int i = 0; i < numberOfPoints.at( 0 ); i++ )
    {
        xDataPoints.push_back( bounds[ 0 ][ 0 ] + static_cast< double >( i ) * xSpacing );
    }

    std::vector< double > yDataPoints;
    for( int j = 0; j < numberOfPoints.at( 1 ); j++ )
    {
        yDataPoints.push_back( bounds[ 0 ][ 1 ] + static_cast< double >( j ) * ySpacing );
    }

    for( int i = 0; i < numberOfPoints.at( 0 ); i++ )
    {
        std::cout<<"Grid search "<<i<<std::endl;
        for( int j = 0; j < numberOfPoints.at( 1 ); j++ )
        {
            decisionVector[ 0 ] = xDataPoints[ i ];
            decisionVector[ 1 ] = yDataPoints[ j ];

            gridSearch( i, j ) = problem.fitness( decisionVector ).at( 0 );
        }
    }

    tudat::input_output::writeMatrixToFile( gridSearch, fileName + ".dat" , 16, getOutputPath( ) );
    tudat::input_output::writeMatrixToFile( tudat::utilities::convertStlVectorToEigenVector(
                                                xDataPoints ), fileName + "_x_data.dat", 16, getOutputPath( ) );
    tudat::input_output::writeMatrixToFile( tudat::utilities::convertStlVectorToEigenVector(
                                                yDataPoints ), fileName + "_y_data.dat", 16, getOutputPath( ) );

}

}

#endif // TUDAT_PAGMO_APPLICATIONOUTPUT_H

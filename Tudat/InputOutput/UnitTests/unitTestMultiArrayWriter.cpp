/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/erase.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/InputOutput/streamFilters.h"
#include "Tudat/Basics/utilities.h"

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/multiDimensionalArrayWriter.h"
#include "Tudat/InputOutput/multiDimensionalArrayReader.h"
#include "Tudat/InputOutput/matrixTextFileReader.h"
#include "Tudat/InputOutput/mapTextFileReader.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_multi_array_writer )

// Test if multi-array file writer is working correctly, if 1 dimension is used
BOOST_AUTO_TEST_CASE( testMultiArrayWriterOneIndependentVariable )
{

    // Test first method of writing 1-dimensional files: one dependent variable per file
    {
        // Define tolerance for equality
        const double equalityTolerance = 1e-8;

        // Read 1-dimensional multi-array
        std::string inFileName1 = tudat::input_output::getTudatRootPath( )
                + "Astrodynamics/Aerodynamics/UnitTests/tabulatedDragCoefficient.txt";
        std::string inFileName2 = tudat::input_output::getTudatRootPath( )
                + "Astrodynamics/Aerodynamics/UnitTests/tabulatedDragCoefficient.txt";

        // Extract drag coefficient data
        Eigen::MatrixXd fileContents1 = tudat::input_output::readMatrixFromFile( inFileName1 );
        Eigen::VectorXd dragCoefficientsBefore = fileContents1.col( 1 );
        Eigen::VectorXd independentVariables1Before = fileContents1.col( 0 );

        // Extract lift coefficient data
        Eigen::MatrixXd fileContents2 = tudat::input_output::readMatrixFromFile( inFileName2 );
        Eigen::VectorXd liftCoefficientsBefore = fileContents2.col( 1 );
        Eigen::VectorXd independentVariables2Before = fileContents2.col( 0 );

        // Check that independent variables are equal
        BOOST_CHECK_EQUAL( independentVariables1Before.size( ), independentVariables2Before.size( ) );
        for ( unsigned int i = 0; i < independentVariables1Before.size( ); i++ )
        {
            BOOST_CHECK_EQUAL( independentVariables1Before[ i ], independentVariables2Before[ i ] );
        }

        // Create multi-array of aerodynamic coefficients
        Eigen::Vector6d currentCoefficients = Eigen::Vector6d::Zero( );
        boost::multi_array< Eigen::Vector6d, 1 > aerodynamicCoefficients( boost::extents[ independentVariables1Before.size( ) ] );
        for ( unsigned int i = 0; i < independentVariables1Before.size( ); i++ )
        {
            currentCoefficients[ 0 ] = dragCoefficientsBefore[ i ];
            currentCoefficients[ 2 ] = liftCoefficientsBefore[ i ];
            aerodynamicCoefficients[ i ] = currentCoefficients;
        }

        // Transform independent variables to vector of vectors
        std::vector< std::vector< double > > independentVariables;
        independentVariables.push_back( tudat::utilities::convertEigenVectorToStlVector< double >( independentVariables1Before ) );

        // Create new multi-array files
        std::map< int, std::string > outFileNamesMap;
        outFileNamesMap[ 0 ] = tudat::input_output::getTudatRootPath( ) + "InputOutput/UnitTests/tabulatedDragCoefficient1_copy.txt";
        outFileNamesMap[ 2 ] = tudat::input_output::getTudatRootPath( ) + "InputOutput/UnitTests/tabulatedDragCoefficient2_copy.txt";

        // Write multi-array to new files
        tudat::input_output::MultiArrayFileWriter< 1, 6 >::writeMultiArrayAndIndependentVariablesToFiles(
                    outFileNamesMap, independentVariables, aerodynamicCoefficients );

        // Read files once more
        // Extract drag coefficient data
        fileContents1 = tudat::input_output::readMatrixFromFile( inFileName1 );
        Eigen::VectorXd dragCoefficientsAfter = fileContents1.col( 1 );
        Eigen::VectorXd independentVariables1After = fileContents1.col( 0 );

        // Extract lift coefficient data
        fileContents2 = tudat::input_output::readMatrixFromFile( inFileName2 );
        Eigen::VectorXd liftCoefficientsAfter = fileContents2.col( 1 );
        Eigen::VectorXd independentVariables2After = fileContents2.col( 0 );

        // Check that independent variables are equal
        BOOST_CHECK_EQUAL( independentVariables1Before.size( ), independentVariables1After.size( ) );
        for ( unsigned int i = 0; i < independentVariables1Before.size( ); i++ )
        {
            BOOST_CHECK_EQUAL( independentVariables1Before[ i ], independentVariables1After[ i ] );
        }
        BOOST_CHECK_EQUAL( independentVariables2Before.size( ), independentVariables2After.size( ) );
        for ( unsigned int i = 0; i < independentVariables1Before.size( ); i++ )
        {
            BOOST_CHECK_EQUAL( independentVariables2Before[ i ], independentVariables2After[ i ] );
        }

        // Check that coefficients are equal
        for ( unsigned int i = 0; i < independentVariables1Before.size( ); i++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( dragCoefficientsAfter[ i ], dragCoefficientsBefore[ i ],
                                        equalityTolerance );
            BOOST_CHECK_CLOSE_FRACTION( liftCoefficientsAfter[ i ], liftCoefficientsBefore[ i ],
                                        equalityTolerance );
        }
    }

    // Test second method of writing 1-dimensional files: independent variable and all dependent variables along columns (side-by-side)
    {
        // Define tolerance for equality
        const double equalityTolerance = 1e-8;

        // Read 1-dimensional multi-array
        std::string inFileName1 = tudat::input_output::getTudatRootPath( )
                + "Astrodynamics/Aerodynamics/UnitTests/tabulatedDragCoefficient.txt";
        std::string inFileName2 = tudat::input_output::getTudatRootPath( )
                + "Astrodynamics/Aerodynamics/UnitTests/tabulatedDragCoefficient.txt";

        // Extract drag coefficient data
        Eigen::MatrixXd fileContents1 = tudat::input_output::readMatrixFromFile( inFileName1 );
        Eigen::VectorXd dragCoefficientsBefore = fileContents1.col( 1 );
        Eigen::VectorXd independentVariables1Before = fileContents1.col( 0 );

        // Extract lift coefficient data
        Eigen::MatrixXd fileContents2 = tudat::input_output::readMatrixFromFile( inFileName2 );
        Eigen::VectorXd liftCoefficientsBefore = fileContents2.col( 1 );
        Eigen::VectorXd independentVariables2Before = fileContents2.col( 0 );

        // Check that independent variables are equal
        BOOST_CHECK_EQUAL( independentVariables1Before.size( ), independentVariables2Before.size( ) );
        for ( unsigned int i = 0; i < independentVariables1Before.size( ); i++ )
        {
            BOOST_CHECK_EQUAL( independentVariables1Before[ i ], independentVariables2Before[ i ] );
        }

        // Create multi-array of aerodynamic coefficients
        Eigen::Vector6d currentCoefficients = Eigen::Vector6d::Zero( );
        boost::multi_array< Eigen::Vector6d, 1 > aerodynamicCoefficients( boost::extents[ independentVariables1Before.size( ) ] );
        for ( unsigned int i = 0; i < independentVariables1Before.size( ); i++ )
        {
            currentCoefficients[ 0 ] = dragCoefficientsBefore[ i ];
            currentCoefficients[ 2 ] = liftCoefficientsBefore[ i ];
            aerodynamicCoefficients[ i ] = currentCoefficients;
        }

        // Transform independent variables to vector of vectors
        std::vector< std::vector< double > > independentVariables;
        independentVariables.push_back( tudat::utilities::convertEigenVectorToStlVector< double >( independentVariables1Before ) );

        // Create new multi-array files
        std::string outFileName = tudat::input_output::getTudatRootPath( ) + "InputOutput/UnitTests/tabulatedDragCoefficient3_copy.txt";

        // Write multi-array to new files
        tudat::input_output::MultiArrayFileWriter< 1, 6 >::writeMultiArrayAndIndependentVariablesToFiles(
                    outFileName, { 0, 2 }, independentVariables, aerodynamicCoefficients );

        // Read files once more
        fileContents1 = tudat::input_output::readMatrixFromFile( outFileName );
        Eigen::VectorXd liftCoefficientsAfter = fileContents1.col( 2 );
        Eigen::VectorXd dragCoefficientsAfter = fileContents1.col( 1 );
        Eigen::VectorXd independentVariables1After = fileContents1.col( 0 );

        // Check that independent variables are equal
        BOOST_CHECK_EQUAL( independentVariables1Before.size( ), independentVariables1After.size( ) );
        for ( unsigned int i = 0; i < independentVariables1Before.size( ); i++ )
        {
            BOOST_CHECK_EQUAL( independentVariables1Before[ i ], independentVariables1After[ i ] );
        }

        // Check that coefficients are equal
        for ( unsigned int i = 0; i < independentVariables1Before.size( ); i++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( dragCoefficientsAfter[ i ], dragCoefficientsBefore[ i ],
                                        equalityTolerance );
            BOOST_CHECK_CLOSE_FRACTION( liftCoefficientsAfter[ i ], liftCoefficientsBefore[ i ],
                                        equalityTolerance );
        }
    }
}

// Test if multi-array file writer is working correctly, when 2 or more dimensions are used
BOOST_AUTO_TEST_CASE( testMultiArrayWriterMultiIndependentVariables )
{
    // Test functionality of 2-dimensional multi-array writer
    {
        // Define tolerance for equality
        const double equalityTolerance = 1e-8;

        // Read 2-dimensional multi-array
        std::string inFileName1 = tudat::input_output::getTudatRootPath( )
                + "Astrodynamics/Aerodynamics/UnitTests/aurora_CD.txt";
        std::string inFileName2 = tudat::input_output::getTudatRootPath( )
                + "Astrodynamics/Aerodynamics/UnitTests/aurora_CL.txt";

        // Extract drag coefficient data
        std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > > fileContents1 =
                tudat::input_output::MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables( inFileName1 );
        boost::multi_array< double, 2 > dragCoefficientsBefore = fileContents1.first;
        std::vector< std::vector< double > > independentVariables1Before = fileContents1.second;

        // Extract lift coefficient data
        std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > > fileContents2 =
                tudat::input_output::MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables( inFileName2 );
        boost::multi_array< double, 2 > liftCoefficientsBefore = fileContents2.first;
        std::vector< std::vector< double > > independentVariables2Before = fileContents2.second;

        // Check that independent variables are equal
        BOOST_CHECK_EQUAL( independentVariables1Before.size( ), independentVariables2Before.size( ) );
        for ( unsigned int i = 0; i < independentVariables1Before.size( ); i++ )
        {
            BOOST_CHECK_EQUAL( independentVariables1Before.at( i ).size( ), independentVariables2Before.at( i ).size( ) );
            for ( unsigned int j = 0; j < independentVariables1Before.at( i ).size( ); j++ )
            {
                BOOST_CHECK_CLOSE_FRACTION( independentVariables1Before.at( i ).at( j ),
                                            independentVariables2Before.at( i ).at( j ),
                                            equalityTolerance );
            }
        }

        // Create multi-array of aerodynamic coefficients
        Eigen::Vector6d currentCoefficients = Eigen::Vector6d::Zero( );
        boost::multi_array< Eigen::Vector6d, 2 > aerodynamicCoefficients(
                    boost::extents[ independentVariables1Before.at( 0 ).size( ) ][ independentVariables1Before.at( 1 ).size( ) ] );
        for ( unsigned int i = 0; i < independentVariables1Before.at( 0 ).size( ); i++ )
        {
            for ( unsigned int j = 0; j < independentVariables1Before.at( 1 ).size( ); j++ )
            {
                currentCoefficients[ 0 ] = dragCoefficientsBefore[ i ][ j ];
                currentCoefficients[ 2 ] = liftCoefficientsBefore[ i ][ j ];
                aerodynamicCoefficients[ i ][ j ] = currentCoefficients;
            }
        }

        // Create new multi-array files
        std::map< int, std::string > outFileNamesMap;
        outFileNamesMap[ 0 ] = tudat::input_output::getTudatRootPath( ) + "InputOutput/UnitTests/aurora_CD_copy.txt";
        outFileNamesMap[ 2 ] = tudat::input_output::getTudatRootPath( ) + "InputOutput/UnitTests/aurora_CL_copy.txt";

        // Write multi-array to new files
        tudat::input_output::MultiArrayFileWriter< 2, 6 >::writeMultiArrayAndIndependentVariablesToFiles(
                    outFileNamesMap, independentVariables1Before, aerodynamicCoefficients );

        // Read files once more
        // Extract drag coefficient data
        fileContents1 = tudat::input_output::MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables( outFileNamesMap[ 0 ] );
        boost::multi_array< double, 2 > dragCoefficientsAfter = fileContents1.first;
        std::vector< std::vector< double > > independentVariables1After = fileContents1.second;

        // Extract lift coefficient data
        fileContents2 = tudat::input_output::MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables( outFileNamesMap[ 2 ] );
        boost::multi_array< double, 2 > liftCoefficientsAfter = fileContents2.first;
        std::vector< std::vector< double > > independentVariables2After = fileContents2.second;

        // Check that independent variables are equal
        BOOST_CHECK_EQUAL( independentVariables1Before.size( ), independentVariables1After.size( ) );
        for ( unsigned int i = 0; i < independentVariables1Before.size( ); i++ )
        {
            BOOST_CHECK_EQUAL( independentVariables1Before.at( i ).size( ), independentVariables1After.at( i ).size( ) );
            for ( unsigned int j = 0; j < independentVariables1Before.at( i ).size( ); j++ )
            {
                BOOST_CHECK_CLOSE_FRACTION( independentVariables1Before.at( i ).at( j ),
                                            independentVariables1After.at( i ).at( j ),
                                            equalityTolerance );
            }
        }
        BOOST_CHECK_EQUAL( independentVariables2Before.size( ), independentVariables2After.size( ) );
        for ( unsigned int i = 0; i < independentVariables1Before.size( ); i++ )
        {
            BOOST_CHECK_EQUAL( independentVariables2Before.at( i ).size( ), independentVariables2After.at( i ).size( ) );
            for ( unsigned int j = 0; j < independentVariables1Before.at( i ).size( ); j++ )
            {
                BOOST_CHECK_CLOSE_FRACTION( independentVariables2Before.at( i ).at( j ),
                                            independentVariables2After.at( i ).at( j ),
                                            equalityTolerance );
            }
        }

        // Check that coefficients are equal
        for ( unsigned int i = 0; i < independentVariables1Before.at( 0 ).size( ); i++ )
        {
            for ( unsigned int j = 0; j < independentVariables1Before.at( 1 ).size( ); j++ )
            {
                BOOST_CHECK_CLOSE_FRACTION( dragCoefficientsAfter[ i ][ j ], dragCoefficientsBefore[ i ][ j ],
                                            equalityTolerance );
                BOOST_CHECK_CLOSE_FRACTION( liftCoefficientsAfter[ i ][ j ], liftCoefficientsBefore[ i ][ j ],
                                            equalityTolerance );
            }
        }
    }

    // Test functionality of 3-dimensional multi-array writer
    {
        // Define tolerance for equality
        const double equalityTolerance = 1e-8;

        // Read 3-dimensional multi-array
        std::string inFileName1 = tudat::input_output::getTudatRootPath( )
                + "Astrodynamics/Aerodynamics/UnitTests/dCDwTest.txt";
        std::string inFileName2 = tudat::input_output::getTudatRootPath( )
                + "Astrodynamics/Aerodynamics/UnitTests/dCDwTest.txt";

        // Extract drag coefficient data
        std::pair< boost::multi_array< double, 3 >, std::vector< std::vector< double > > > fileContents1 =
                tudat::input_output::MultiArrayFileReader< 3 >::readMultiArrayAndIndependentVariables( inFileName1 );
        boost::multi_array< double, 3 > dragCoefficientsBefore = fileContents1.first;
        std::vector< std::vector< double > > independentVariables1Before = fileContents1.second;

        // Extract lift coefficient data
        std::pair< boost::multi_array< double, 3 >, std::vector< std::vector< double > > > fileContents2 =
                tudat::input_output::MultiArrayFileReader< 3 >::readMultiArrayAndIndependentVariables( inFileName2 );
        boost::multi_array< double, 3 > liftCoefficientsBefore = fileContents2.first;
        std::vector< std::vector< double > > independentVariables2Before = fileContents2.second;

        // Check that independent variables are equal
        BOOST_CHECK_EQUAL( independentVariables1Before.size( ), independentVariables2Before.size( ) );
        for ( unsigned int i = 0; i < independentVariables1Before.size( ); i++ )
        {
            BOOST_CHECK_EQUAL( independentVariables1Before.at( i ).size( ), independentVariables2Before.at( i ).size( ) );
            for ( unsigned int j = 0; j < independentVariables1Before.at( i ).size( ); j++ )
            {
                BOOST_CHECK_CLOSE_FRACTION( independentVariables1Before.at( i ).at( j ),
                                            independentVariables2Before.at( i ).at( j ),
                                            equalityTolerance );
            }
        }

        // Create multi-array of aerodynamic coefficients
        Eigen::Vector6d currentCoefficients = Eigen::Vector6d::Zero( );
        boost::multi_array< Eigen::Vector6d, 3 > aerodynamicCoefficients(
                    boost::extents[ independentVariables1Before.at( 0 ).size( ) ][
                independentVariables1Before.at( 1 ).size( ) ][ independentVariables1Before.at( 2 ).size( ) ] );
        for ( unsigned int i = 0; i < independentVariables1Before.at( 0 ).size( ); i++ )
        {
            for ( unsigned int j = 0; j < independentVariables1Before.at( 1 ).size( ); j++ )
            {
                for ( unsigned int k = 0; k < independentVariables1Before.at( 2 ).size( ); k++ )
                {
                    currentCoefficients[ 0 ] = dragCoefficientsBefore[ i ][ j ][ k ];
                    currentCoefficients[ 2 ] = liftCoefficientsBefore[ i ][ j ][ k ];
                    aerodynamicCoefficients[ i ][ j ][ k ] = currentCoefficients;
                }
            }
        }

        // Create new multi-array files
        std::map< int, std::string > outFileNamesMap;
        outFileNamesMap[ 0 ] = tudat::input_output::getTudatRootPath( ) + "InputOutput/UnitTests/dCDwTest1_copy.txt";
        outFileNamesMap[ 2 ] = tudat::input_output::getTudatRootPath( ) + "InputOutput/UnitTests/dCDwTest2_copy.txt";

        // Write multi-array to new files
        tudat::input_output::MultiArrayFileWriter< 3, 6 >::writeMultiArrayAndIndependentVariablesToFiles(
                    outFileNamesMap, independentVariables1Before, aerodynamicCoefficients );

        // Read files once more
        // Extract drag coefficient data
        fileContents1 = tudat::input_output::MultiArrayFileReader< 3 >::readMultiArrayAndIndependentVariables( outFileNamesMap[ 0 ] );
        boost::multi_array< double, 3 > dragCoefficientsAfter = fileContents1.first;
        std::vector< std::vector< double > > independentVariables1After = fileContents1.second;

        // Extract lift coefficient data
        fileContents2 = tudat::input_output::MultiArrayFileReader< 3 >::readMultiArrayAndIndependentVariables( outFileNamesMap[ 2 ] );
        boost::multi_array< double, 3 > liftCoefficientsAfter = fileContents2.first;
        std::vector< std::vector< double > > independentVariables2After = fileContents2.second;

        // Check that independent variables are equal
        BOOST_CHECK_EQUAL( independentVariables1Before.size( ), independentVariables1After.size( ) );
        for ( unsigned int i = 0; i < independentVariables1Before.size( ); i++ )
        {
            BOOST_CHECK_EQUAL( independentVariables1Before.at( i ).size( ), independentVariables1After.at( i ).size( ) );
            for ( unsigned int j = 0; j < independentVariables1Before.at( i ).size( ); j++ )
            {
                BOOST_CHECK_CLOSE_FRACTION( independentVariables1Before.at( i ).at( j ),
                                            independentVariables1After.at( i ).at( j ),
                                            equalityTolerance );
            }
        }
        BOOST_CHECK_EQUAL( independentVariables2Before.size( ), independentVariables2After.size( ) );
        for ( unsigned int i = 0; i < independentVariables1Before.size( ); i++ )
        {
            BOOST_CHECK_EQUAL( independentVariables2Before.at( i ).size( ), independentVariables2After.at( i ).size( ) );
            for ( unsigned int j = 0; j < independentVariables1Before.at( i ).size( ); j++ )
            {
                BOOST_CHECK_CLOSE_FRACTION( independentVariables2Before.at( i ).at( j ),
                                            independentVariables2After.at( i ).at( j ),
                                            equalityTolerance );
            }
        }

        // Check that coefficients are equal
        for ( unsigned int i = 0; i < independentVariables1Before.at( 0 ).size( ); i++ )
        {
            for ( unsigned int j = 0; j < independentVariables1Before.at( 1 ).size( ); j++ )
            {
                for ( unsigned int k = 0; k < independentVariables1Before.at( 2 ).size( ); k++ )
                {
                    BOOST_CHECK_CLOSE_FRACTION( dragCoefficientsAfter[ i ][ j ][ k ], dragCoefficientsBefore[ i ][ j ][ k ],
                                                equalityTolerance );
                    BOOST_CHECK_CLOSE_FRACTION( liftCoefficientsAfter[ i ][ j ][ k ], liftCoefficientsBefore[ i ][ j ][ k ],
                                                equalityTolerance );
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat


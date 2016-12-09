#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

#include "Tudat/InputOutput/multiDimensionalArrayReader.h"

namespace tudat
{

namespace input_output
{

boost::multi_array< double, 1 > readOneDimensionalCoefficientFile(
        const std::vector< int > independentVariableSize,
        const Eigen::MatrixXd& coefficientsBlock )
{
    boost::multi_array< double, 1 > coefficientMultiarray;

    if( independentVariableSize.size( ) == 1 )
    {
        coefficientMultiarray.resize( boost::extents[ independentVariableSize.at( 0 ) ] );

        for( int i = 0; i < independentVariableSize.at( 0 ); i++ )
        {
            coefficientMultiarray[ i ] = coefficientsBlock( i, 0 );
        }
    }
    return coefficientMultiarray;
}

boost::multi_array< double, 2 > readTwoDimensionalCoefficientFile(
        const std::vector< int > independentVariableSize,
        const Eigen::MatrixXd& coefficientsBlock )
{
    boost::multi_array< double, 2 > coefficientMultiarray;

    if( independentVariableSize.size( ) == 2 )
    {
        coefficientMultiarray.resize( boost::extents[ independentVariableSize.at( 0 ) ][ independentVariableSize.at( 1 ) ]  );

        for( int i = 0; i < independentVariableSize.at( 0 ); i++ )
        {
            for( int j = 0; j < independentVariableSize.at( 1 ); j++ )
            {
                coefficientMultiarray[ i ][ j ] = coefficientsBlock( i, j );
            }
        }
    }
    return coefficientMultiarray;
}

boost::multi_array< double, 3 > readThreeDimensionalCoefficientFile(
        const std::vector< int > independentVariableSize,
        const Eigen::MatrixXd& coefficientsBlock )
{
    boost::multi_array< double, 3 > coefficientMultiarray;

    if( independentVariableSize.size( ) == 3 )
    {
        coefficientMultiarray.resize( boost::extents[ independentVariableSize.at( 0 ) ][ independentVariableSize.at( 1 ) ][ independentVariableSize.at( 2 ) ] );
        int currentStartRow = 0;
        for( int k = 0; k < independentVariableSize.at( 2 ); k++ )
        {
            for( int i = 0; i < independentVariableSize.at( 0 ); i++ )
            {
                for( int j = 0; j < independentVariableSize.at( 1 ); j++ )
                {
                    coefficientMultiarray[ i ][ j ][ k ] = coefficientsBlock( i + currentStartRow, j );
                }
            }
            currentStartRow += independentVariableSize.at( 0 );
        }
    }
    return coefficientMultiarray;
}

void readCoefficientsFile(
        const std::string fileName,
        std::vector< std::vector< double > >& independentVariables,
        Eigen::MatrixXd& coefficientBlock )
{
    // Open file and create file stream.
    std::fstream stream( fileName.c_str( ), std::ios::in );

    // Check if file opened correctly.
    if ( stream.fail( ) )
    {
        boost::throw_exception( std::runtime_error( boost::str(
                                                        boost::format( "Data file '%s' could not be opened." ) % fileName.c_str( ) ) ) );
    }

    // Initialize boolean that gets set to true once the file header is passed.
    bool isHeaderPassed = 0;
    bool isFirstLinePassed = 0;


    // Line based parsing
    std::string line;
    std::vector< std::string > vectorOfIndividualStrings;

    int numberOfDataLinesParsed = 0;
    int numberOfIndependentVariables = -1;
    while ( !stream.fail( ) && !stream.eof( ) )
    {
        // Get line from stream
        std::getline( stream, line );

        // Trim input string (removes all leading and trailing whitespaces).
        boost::algorithm::trim( line );

        if( line.size( ) > 0 && !( line.at( 0 ) == '#' ) )
        {
            // Split string into multiple strings, each containing one element from a line from the
            // data file.
            boost::algorithm::split( vectorOfIndividualStrings,
                                     line,
                                     boost::algorithm::is_any_of( "\t ;, " ),
                                     boost::algorithm::token_compress_on );
            if( !isFirstLinePassed )
            {
                if( vectorOfIndividualStrings.size( ) != 1 )
                {
                    throw std::runtime_error( "Error when reading multi-array, expected number of independent variables" );
                }
                numberOfIndependentVariables = boost::lexical_cast< int >( vectorOfIndividualStrings.at( 0 ) );
                isFirstLinePassed = true;
            }
            else if( !isHeaderPassed )
            {
                std::vector< double > currentDataPoints;
                for( unsigned int i = 0; i < vectorOfIndividualStrings.size( ); i++ )
                {
                    currentDataPoints.push_back( boost::lexical_cast< double >( vectorOfIndividualStrings.at( i ) ) );
                }
                independentVariables.push_back( currentDataPoints );
            }
            else if( isHeaderPassed )
            {
                if( vectorOfIndividualStrings.size( ) != static_cast< unsigned int >( coefficientBlock.cols( ) ) )
                {
                    throw std::runtime_error( "Error on data line " + boost::lexical_cast< std::string >( numberOfDataLinesParsed ) +
                                              " found " + boost::lexical_cast< std::string >( vectorOfIndividualStrings.size( ) ) +
                                              " columns, but expected " +  boost::lexical_cast< std::string >( coefficientBlock.cols( ) ) );
                }
                else if( numberOfDataLinesParsed > coefficientBlock.rows( ) )
                {
                    throw std::runtime_error( "Error on data line " + boost::lexical_cast< std::string >( numberOfDataLinesParsed ) +
                                              " expected " +  boost::lexical_cast< std::string >( coefficientBlock.rows( ) ) + "rows" );
                }
                else
                {
                    for( unsigned int i = 0; i < vectorOfIndividualStrings.size( ); i++ )
                    {
                        coefficientBlock( numberOfDataLinesParsed, i ) = boost::lexical_cast< double >( vectorOfIndividualStrings.at( i ) );
                    }
                    numberOfDataLinesParsed++;
                }
            }

            if( ( static_cast< int >( independentVariables.size( ) ) == numberOfIndependentVariables ) && !isHeaderPassed )
            {
                if( independentVariables.size( ) == 0 )
                {
                    throw std::runtime_error( "Error when reading multi-array, no header found" );
                }
                else
                {
                    int numberOfRows = independentVariables.at( 0 ).size( );
                    int numberOfColumns = 1;
                    if( independentVariables.size( ) > 1 )
                    {
                        numberOfColumns = independentVariables.at( 1 ).size( );
                        for( unsigned int i = 2; i < independentVariables.size( ); i++ )
                        {
                            numberOfRows *= independentVariables.at( i ).size( );
                        }
                    }
                    coefficientBlock.setZero( numberOfRows, numberOfColumns );
                }
                isHeaderPassed = true;
            }

        }
    }

    if( numberOfDataLinesParsed != coefficientBlock.rows( ) )
    {
        throw std::runtime_error( "Error at end of coefficient file reader, found " + boost::lexical_cast< std::string >( numberOfDataLinesParsed ) +
                                  " lines, but expected " +  boost::lexical_cast< std::string >( coefficientBlock.rows( ) ) + "rows" );
    }
}

}

}

/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_READHISTORYFROMFILE_H
#define TUDAT_READHISTORYFROMFILE_H

#include <map>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <Eigen/Core>

namespace tudat
{

namespace input_output
{

//! Function to read a time history of Eigen MatrixXd data from a file
/*!
 *  Function to read a time history of Eigen MatrixXd data from a file, as a map with time (key) and associated MatrixXd (value)
 *  \param singleMatrixRows Number of rows per matrix
 *  \param singleMatrixColumns Number of columns per matrix
 *  \param fileName File name to load
 *  \return Matrix history from file.
 */
template< typename TimeType, typename StateScalarType >
std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > readMatrixHistoryFromFile(
        const int singleMatrixRows, const int singleMatrixColumns, const std::string& fileName )
{
    std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > matrixHistory;

    std::ifstream fileStream;
    fileStream.open( fileName );
    if ( !fileStream.is_open( ) )
    {
        throw std::runtime_error( "Data file: " + fileName + " could not be opened." );
    }

    TimeType currentTime;
    double doubleCurrentTime;
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > currentMatrix =
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >::Zero( singleMatrixRows, singleMatrixColumns );

    while( !fileStream.eof( ) )
    {
        fileStream >> doubleCurrentTime;
        currentTime = doubleCurrentTime;
        for( int i = 0; i < singleMatrixRows; i++ )
        {
            for( int j = 0; j < singleMatrixColumns; j++ )
            {
                fileStream >> currentMatrix( i, j );
            }
        }

        matrixHistory[ currentTime ] = currentMatrix;

    }
    fileStream.close( );

    return matrixHistory;
}

//! Function to read a time history of Eigen VectorXd data from a file
/*!
 *  Function to read a time history of Eigen VectorXd data from a file, as a map with time (key) and associated VectorXd (value)
 *  \param singleMatrixRows Number of rows per vector
 *  \param fileName File name to load
 *  \return Vector history from file.
 */
template< typename TimeType, typename StateScalarType >
std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > readVectorHistoryFromFile(
        const int singleMatrixRows, const std::string& fileName )
{
    std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > matrixHistory;

    std::ifstream fileStream;
    fileStream.open( fileName );
    if ( !fileStream.is_open( ) )
    {
        throw std::runtime_error( "Data file: " + fileName + " could not be opened." );
    }

    TimeType currentTime;
    double doubleCurrentTime;
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentMatrix =
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( singleMatrixRows, 1 );

    while( !fileStream.eof( ) )
    {
        fileStream >> doubleCurrentTime;
        currentTime = doubleCurrentTime;
        for( int i = 0; i < singleMatrixRows; i++ )
        {
            fileStream >> currentMatrix( i, 0 );
        }

        matrixHistory[ currentTime ] = currentMatrix;

    }
    fileStream.close( );

    return matrixHistory;
}

//! Function to read a time history of scalar data from a file
/*!
 *  Function to read a time history of scalar data from a file, as a map with time (key) and associated scalar (value)
 *  \param fileName File name to load
 *  \return Scalar history from file.
 */
template< typename S, typename T >
std::map< S, T > readScalarHistoryFromFile( const std::string& fileName )
{
    std::map< S, T > dataMap;
    std::ifstream fileStream;
    fileStream.open( fileName );
    S key;
    T value;
    while( !fileStream.eof( ) )
    {
        fileStream >> key;
        fileStream >> value;

        dataMap[ key ] = value;
    }
    fileStream.close( );
    return dataMap;
}


}

}

#endif // TUDAT_READHISTORYFROMFILE_H

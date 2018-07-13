#ifndef READHISTORYFROMFILE_H
#define READHISTORYFROMFILE_H

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

template< typename TimeType, typename StateScalarType >
std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > readMatrixHistoryFromFile(
        const int singleMatrixRows, const int singleMatrixColumns, const std::string& fileName )
{
    std::cout<<"Reading file"<<std::endl;
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
    std::cout<<"File read"<<std::endl;

    return matrixHistory;
}

template< typename TimeType, typename StateScalarType >
std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > readVectorHistoryFromFile(
        const int singleMatrixRows, const std::string& fileName )
{
    std::cout<<"Reading file"<<std::endl;
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
    std::cout<<"File read"<<std::endl;

    return matrixHistory;
}

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

#endif // READHISTORYFROMFILE_H

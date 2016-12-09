#ifndef MULTIDIMENSIONALARRAYREADER_H
#define MULTIDIMENSIONALARRAYREADER_H

#include <boost/multi_array.hpp>

#include <Eigen/Core>

namespace tudat
{

namespace input_output
{

boost::multi_array< double, 1 > readOneDimensionalCoefficientFile(
        const std::vector< int > independentVariableSize,
        const Eigen::MatrixXd& coefficientsBlock );

boost::multi_array< double, 2 > readTwoDimensionalCoefficientFile(
        const std::vector< int > independentVariableSize,
        const Eigen::MatrixXd& coefficientsBlock );

boost::multi_array< double, 3 > readThreeDimensionalCoefficientFile(
        const std::vector< int > independentVariableSize,
        const Eigen::MatrixXd& coefficientsBlock );

void readCoefficientsFile(
        const std::string fileName,
        std::vector< std::vector< double > >& independentVariables,
        Eigen::MatrixXd& coefficientBlock );

template< int NumberOfDimensions >
class MultiArrayFileReader
{
public:


    static boost::multi_array< double, NumberOfDimensions > readMultiArray(
            const std::string fileName );

    static std::pair< boost::multi_array< double, 1 >, std::vector< std::vector< double > > >
        readMultiArrayAndIndependentVariables(
            const std::string fileName );
};

template< >
class MultiArrayFileReader< 1 >
{
public:

    static boost::multi_array< double, 1 > readMultiArray(
            const std::string fileName )
    {
        return readMultiArrayAndIndependentVariables( fileName ).first;
    }

    static std::pair< boost::multi_array< double, 1 >, std::vector< std::vector< double > > >
        readMultiArrayAndIndependentVariables(
            const std::string fileName )
    {
        std::vector< std::vector< double > > independentVariables;
        Eigen::MatrixXd coefficientBlock;
        readCoefficientsFile( fileName, independentVariables, coefficientBlock );

        if( !( independentVariables.size( ) == 1 ) )
        {
            throw std::runtime_error( "Error when reading 3-D multi-array, wrong number of independent variables found" );
        }

        std::vector< int > independentVariableSize;
        independentVariableSize.push_back( independentVariables.at( 0 ).size( ) );

        return std::make_pair( readOneDimensionalCoefficientFile( independentVariableSize, coefficientBlock ), independentVariables );
    }
};

template< >
class MultiArrayFileReader< 2 >
{
public:

    static boost::multi_array< double, 2 > readMultiArray(
            const std::string fileName )
    {
        return readMultiArrayAndIndependentVariables( fileName ).first;
    }

    static std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > >
        readMultiArrayAndIndependentVariables(
            const std::string fileName )
    {
        std::vector< std::vector< double > > independentVariables;
        Eigen::MatrixXd coefficientBlock;
        readCoefficientsFile( fileName, independentVariables, coefficientBlock );

        if( !( independentVariables.size( ) == 2 ) )
        {
            throw std::runtime_error( "Error when reading 3-D multi-array, wrong number of independent variables found" );
        }

        std::vector< int > independentVariableSize;
        independentVariableSize.push_back( independentVariables.at( 0 ).size( ) );
        independentVariableSize.push_back( independentVariables.at( 1 ).size( ) );

        return std::make_pair( readTwoDimensionalCoefficientFile( independentVariableSize, coefficientBlock ), independentVariables );
    }
};

template< >
class MultiArrayFileReader< 3 >
{
public:

    static boost::multi_array< double, 3 > readMultiArray(
            const std::string fileName )
    {
        return readMultiArrayAndIndependentVariables( fileName ).first;
    }

    static std::pair< boost::multi_array< double, 3 >, std::vector< std::vector< double > > >
        readMultiArrayAndIndependentVariables(
            const std::string fileName )
    {
        std::vector< std::vector< double > > independentVariables;
        Eigen::MatrixXd coefficientBlock;
        readCoefficientsFile( fileName, independentVariables, coefficientBlock );

        if( !( independentVariables.size( ) == 3 ) )
        {
            throw std::runtime_error( "Error when reading 3-D multi-array, wrong number of independent variables found" );
        }

        std::vector< int > independentVariableSize;
        independentVariableSize.push_back( independentVariables.at( 0 ).size( ) );
        independentVariableSize.push_back( independentVariables.at( 1 ).size( ) );
        independentVariableSize.push_back( independentVariables.at( 2 ).size( ) );

        return std::make_pair( readThreeDimensionalCoefficientFile( independentVariableSize, coefficientBlock ), independentVariables );
    }
    };

}

}

#endif // MULTIDIMENSIONALARRAYREADER_H

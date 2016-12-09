#include <map>

#include "Tudat/InputOutput/multiDimensionalArrayReader.h"

namespace tudat
{

namespace input_output
{

boost::multi_array< Eigen::Vector3d, 1 > mergeOneDimensionalCoefficients(
        const boost::multi_array< double, 1 > xComponents,
        const boost::multi_array< double, 1 > yComponents,
        const boost::multi_array< double, 1 > zComponents );

boost::multi_array< Eigen::Vector3d, 2 > mergeTwoDimensionalCoefficients(
        const boost::multi_array< double, 2 > xComponents,
        const boost::multi_array< double, 2 > yComponents,
        const boost::multi_array< double, 2 > zComponents );

boost::multi_array< Eigen::Vector3d, 3 > mergeThreeDimensionalCoefficients(
        const boost::multi_array< double, 3 > xComponents,
        const boost::multi_array< double, 3 > yComponents,
        const boost::multi_array< double, 3 > zComponents );

bool compareIndependentVariables( const std::vector< std::vector< double > >& list1,
                                  const std::vector< std::vector< double > >& list2 );

template< int NumberOfDimensions >
class AerodynamicCoefficientReader
{
public:


    static boost::multi_array< double, NumberOfDimensions > readMultiArray(
            const std::string fileName );

    static std::pair< boost::multi_array< double, NumberOfDimensions >, std::vector< std::vector< double > > >
        readMultiArrayAndIndependentVariables(
            const std::string fileName );
};

template< >
class AerodynamicCoefficientReader< 1 >
{
public:

    static std::pair< boost::multi_array< Eigen::Vector3d, 1 >, std::vector< std::vector< double > > >
        readAerodynamicCoefficients( const std::vector< std::string >& fileNames )
    {
        if( fileNames.size( ) != 3 )
        {
            throw std::runtime_error( "Error when reading 1-Dimensional aeroynamic coefficients, wrong number of files" );
        }

        std::vector< boost::multi_array< double, 1 > > coefficientArrays;
        std::vector< std::vector< double > > independentVariables;
        for( unsigned int i = 0; i < 3; i++ )
        {
            std::pair< boost::multi_array< double, 1 >, std::vector< std::vector< double > > > currentCoefficients =
                    MultiArrayFileReader< 1 >::readMultiArrayAndIndependentVariables( fileNames.at( i ) );
            if( i == 0 )
            {
                independentVariables = currentCoefficients.second;
            }
            else
            {
                bool areIndependentVariablesEqual = compareIndependentVariables(
                            independentVariables, currentCoefficients.second );
                if( !areIndependentVariablesEqual )
                {
                    throw std::runtime_error( "Error when reading 1-Dimensional aeroynamic coefficients, inconsistent aerodynamic coefficients" );
                }
            }
            coefficientArrays.push_back( currentCoefficients.first );
        }
        return std::make_pair(
                    mergeOneDimensionalCoefficients(
                        coefficientArrays.at( 0 ), coefficientArrays.at( 1 ), coefficientArrays.at( 2 ) ), independentVariables );
    }
};

template< >
class AerodynamicCoefficientReader< 2 >
{
public:

    static std::pair< boost::multi_array< Eigen::Vector3d, 2 >, std::vector< std::vector< double > > >
        readAerodynamicCoefficients( const std::vector< std::string >& fileNames )
    {
        if( fileNames.size( ) != 3 )
        {
            throw std::runtime_error( "Error when reading 1-Dimensional aeroynamic coefficients, wrong number of files" );
        }

        std::vector< boost::multi_array< double, 2 > > coefficientArrays;
        std::vector< std::vector< double > > independentVariables;
        for( unsigned int i = 0; i < 3; i++ )
        {
            std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > > currentCoefficients =
                    MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables( fileNames.at( i ) );

            if( i == 0 )
            {
                independentVariables = currentCoefficients.second;
            }
            else
            {
                bool areIndependentVariablesEqual = compareIndependentVariables(
                            independentVariables, currentCoefficients.second );
                if( !areIndependentVariablesEqual )
                {
                    throw std::runtime_error( "Error when reading 1-Dimensional aeroynamic coefficients, inconsistent aerodynamic coefficients" );
                }
            }
            coefficientArrays.push_back( currentCoefficients.first );
        }
        return std::make_pair(
                    mergeTwoDimensionalCoefficients(
                        coefficientArrays.at( 0 ), coefficientArrays.at( 1 ), coefficientArrays.at( 2 ) ), independentVariables );
    }

};

template< >
class AerodynamicCoefficientReader< 3 >
{
public:

    static std::pair< boost::multi_array< Eigen::Vector3d, 3 >, std::vector< std::vector< double > > >
        readAerodynamicCoefficients( const std::vector< std::string >& fileNames )
    {
        if( fileNames.size( ) != 3 )
        {
            throw std::runtime_error( "Error when reading 1-Dimensional aeroynamic coefficients, wrong number of files" );
        }

        std::vector< boost::multi_array< double, 3 > > coefficientArrays;
        std::vector< std::vector< double > > independentVariables;
        for( unsigned int i = 0; i < 3; i++ )
        {
            std::pair< boost::multi_array< double, 3 >, std::vector< std::vector< double > > > currentCoefficients =
                    MultiArrayFileReader< 3 >::readMultiArrayAndIndependentVariables( fileNames.at( i ) );
            if( i == 0 )
            {
                independentVariables = currentCoefficients.second;
            }
            else
            {
                bool areIndependentVariablesEqual = compareIndependentVariables(
                            independentVariables, currentCoefficients.second );
                if( !areIndependentVariablesEqual )
                {
                    throw std::runtime_error( "Error when reading 1-Dimensional aeroynamic coefficients, inconsistent aerodynamic coefficients" );
                }
            }
            coefficientArrays.push_back( currentCoefficients.first );
        }
        return std::make_pair(
                    mergeThreeDimensionalCoefficients(
                        coefficientArrays.at( 0 ), coefficientArrays.at( 1 ), coefficientArrays.at( 2 ) ), independentVariables );
    }

};

}

}

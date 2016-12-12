/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <map>

#include "Tudat/InputOutput/multiDimensionalArrayReader.h"

namespace tudat
{

namespace input_output
{

//! Function to merge three double multi-arrays of 1 dimension into a single Vector3d multi-array
/*!
 *  Function to merge three double multi-arrays of 1 dimension into a single Vector3d multi-array, where the three
 *  double multi-arrays represent the x-, y- and z-components of the Vector3ds.
 *  \param xComponents Multi-array containing the x-components of the Vector3d
 *  \param yComponents Multi-array containing the y-components of the Vector3d
 *  \param zComponents Multi-array containing the z-components of the Vector3d
 *  \return Single multi-array containing Vector3ds according to double multi-arrays.
 */
boost::multi_array< Eigen::Vector3d, 1 > mergeOneDimensionalCoefficients(
        const boost::multi_array< double, 1 > xComponents,
        const boost::multi_array< double, 1 > yComponents,
        const boost::multi_array< double, 1 > zComponents );

//! Function to merge three double multi-arrays of 2 dimension into a single Vector3d multi-array
/*!
 *  Function to merge three double multi-arrays of 2 dimension into a single Vector3d multi-array, where the three
 *  double multi-arrays represent the x-, y- and z-components of the Vector3ds.
 *  \param xComponents Multi-array containing the x-components of the Vector3d
 *  \param yComponents Multi-array containing the y-components of the Vector3d
 *  \param zComponents Multi-array containing the z-components of the Vector3d
 *  \return Single multi-array containing Vector3ds according to double multi-arrays.
 */
boost::multi_array< Eigen::Vector3d, 2 > mergeTwoDimensionalCoefficients(
        const boost::multi_array< double, 2 > xComponents,
        const boost::multi_array< double, 2 > yComponents,
        const boost::multi_array< double, 2 > zComponents );

//! Function to merge three double multi-arrays of 3 dimension into a single Vector3d multi-array
/*!
 *  Function to merge three double multi-arrays of 3 dimension into a single Vector3d multi-array, where the three
 *  double multi-arrays represent the x-, y- and z-components of the Vector3ds.
 *  \param xComponents Multi-array containing the x-components of the Vector3d
 *  \param yComponents Multi-array containing the y-components of the Vector3d
 *  \param zComponents Multi-array containing the z-components of the Vector3d
 *  \return Single multi-array containing Vector3ds according to double multi-arrays.
 */
boost::multi_array< Eigen::Vector3d, 3 > mergeThreeDimensionalCoefficients(
        const boost::multi_array< double, 3 > xComponents,
        const boost::multi_array< double, 3 > yComponents,
        const boost::multi_array< double, 3 > zComponents );

//! Function to compare if two lists of aerodynamic coefficient independent variables are equal
/*!
 * Function to compare if two lists of aerodynamic coefficient independent variables (vector of vector of doubles) are equal
 * \param list1 First list that is to be compared.
 * \param list2 Second list that is to be compared.
 * \return True of the two lists are completely equal in size and contents, false otherwise.
 */
bool compareIndependentVariables( const std::vector< std::vector< double > >& list1,
                                  const std::vector< std::vector< double > >& list2 );

//! Interface class for reading aerodynamic coefficients from files
/*!
 *  Interface class for reading aerodynamic coefficients from files.
 *  This class is used instead of a single templated free function to allow multi-arrays of different sizes to be created
 *  using the same interface. NOTE: The possibility of using a single templated implementation for arbitrary multi-array
 * size should be investigated in the future.
 */
template< int NumberOfDimensions >
class AerodynamicCoefficientReader
{
public:

    //! Function to read a list of aerodynamic coefficients and associated independent variables from a set of files
    /*!
     *  Function to read a list of aerodynamic coefficients and associated independent variables from a set of files.
     *  \param fileNames Vector of size 3, containing the file names for the x-, y- and z- components of the aerodynamic
     *  coefficients. Note that the independent variables for each components must be identical.
     *  \return  Pair: first entry containing multi-array of aerodynamic coefficients, second containing list of independent
     *  variables at which coefficients are defined.
     */
    static std::pair< boost::multi_array< Eigen::Vector3d, NumberOfDimensions >, std::vector< std::vector< double > > >
        readAerodynamicCoefficients( const std::vector< std::string >& fileNames );
};

//! Interface class for reading aerodynamic coefficients of 1 independent variables from files
template< >
class AerodynamicCoefficientReader< 1 >
{
public:

    //! Function to read a list of aerodynamic coefficients and associated independent variables from a set of files
    /*!
     *  Function to read a list of aerodynamic coefficients of 1 independent variables and associated independent variables
     *  from a set of files.
     *  \param fileNames Vector of size 3, containing the file names for the x-, y- and z- components of the aerodynamic
     *  coefficients. Note that the independent variables for each components must be identical.
     *  \return  Pair: first entry containing multi-array of aerodynamic coefficients, second containing list of independent
     *  variables at which coefficients are defined.
     */
    static std::pair< boost::multi_array< Eigen::Vector3d, 1 >, std::vector< std::vector< double > > >
        readAerodynamicCoefficients( const std::vector< std::string >& fileNames )
    {
        if( fileNames.size( ) != 3 )
        {
            throw std::runtime_error( "Error when reading 1-Dimensional aeroynamic coefficients, wrong number of files" );
        }

        // Read data from files.
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

        // Merge coefficient entries
        return std::make_pair(
                    mergeOneDimensionalCoefficients(
                        coefficientArrays.at( 0 ), coefficientArrays.at( 1 ), coefficientArrays.at( 2 ) ), independentVariables );
    }
};

//! Interface class for reading aerodynamic coefficients of 2 independent variables from files
template< >
class AerodynamicCoefficientReader< 2 >
{
public:

    //! Function to read a list of aerodynamic coefficients and associated independent variables from a set of files
    /*!
     *  Function to read a list of aerodynamic coefficients of 2 independent variables and associated independent variables
     *  from a set of files.
     *  \param fileNames Vector of size 3, containing the file names for the x-, y- and z- components of the aerodynamic
     *  coefficients. Note that the independent variables for each components must be identical.
     *  \return  Pair: first entry containing multi-array of aerodynamic coefficients, second containing list of independent
     *  variables at which coefficients are defined.
     */
    static std::pair< boost::multi_array< Eigen::Vector3d, 2 >, std::vector< std::vector< double > > >
        readAerodynamicCoefficients( const std::vector< std::string >& fileNames )
    {
        if( fileNames.size( ) != 3 )
        {
            throw std::runtime_error( "Error when reading 1-Dimensional aeroynamic coefficients, wrong number of files" );
        }

        // Read data from files.
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

        // Merge coefficient entries
        return std::make_pair(
                    mergeTwoDimensionalCoefficients(
                        coefficientArrays.at( 0 ), coefficientArrays.at( 1 ), coefficientArrays.at( 2 ) ), independentVariables );
    }

};

//! Interface class for reading aerodynamic coefficients of 3 independent variables from files
template< >
class AerodynamicCoefficientReader< 3 >
{
public:

    //! Function to read a list of aerodynamic coefficients and associated independent variables from a set of files
    /*!
     *  Function to read a list of aerodynamic coefficients of 3 independent variables and associated independent variables
     *  from a set of files.
     *  \param fileNames Vector of size 3, containing the file names for the x-, y- and z- components of the aerodynamic
     *  coefficients. Note that the independent variables for each components must be identical.
     *  \return  Pair: first entry containing multi-array of aerodynamic coefficients, second containing list of independent
     *  variables at which coefficients are defined.
     */
    static std::pair< boost::multi_array< Eigen::Vector3d, 3 >, std::vector< std::vector< double > > >
        readAerodynamicCoefficients( const std::vector< std::string >& fileNames )
    {
        if( fileNames.size( ) != 3 )
        {
            throw std::runtime_error( "Error when reading 1-Dimensional aeroynamic coefficients, wrong number of files" );
        }

        // Read data from files.
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

        // Merge coefficient entries
        return std::make_pair(
                    mergeThreeDimensionalCoefficients(
                        coefficientArrays.at( 0 ), coefficientArrays.at( 1 ), coefficientArrays.at( 2 ) ), independentVariables );
    }

};

}

}

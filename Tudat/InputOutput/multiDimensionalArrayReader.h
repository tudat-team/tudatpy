/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MULTIDIMENSIONALARRAYREADER_H
#define TUDAT_MULTIDIMENSIONALARRAYREADER_H

#include <boost/multi_array.hpp>

#include <Eigen/Core>

namespace tudat
{

namespace input_output
{

//! Function to parse a block of values read from a file into a multi-array of size 1.
/*!
 * Function to parse a block of values read from a file into a multi-array of size 1.
 * The size of the independent variable is provided as input from this file.
 * \param independentVariableSize Size of independent variables (must be of size 1)
 * \param coefficientsBlock Raw number block as read from file.
 * \return Multi-array of size 1, according to input data.
 */
boost::multi_array< double, 1 > parseRawOneDimensionalCoefficientsFromFile(
        const std::vector< int > independentVariableSize,
        const Eigen::MatrixXd& coefficientsBlock );

//! Function to parse a block of values read from a file into a multi-array of size 2.
/*!
 * Function to parse a block of values read from a file into a multi-array of size 2.
 * The size of the independent variables is provided as input from this file.
 * \param independentVariableSize Size of independent variables (must be of size 2)
 * \param coefficientsBlock Raw number block as read from file.
 * \return Multi-array of size 2, according to input data.
 */
boost::multi_array< double, 2 > parseRawTwoDimensionalCoefficientsFromFile(
        const std::vector< int > independentVariableSize,
        const Eigen::MatrixXd& coefficientsBlock );

//! Function to parse a block of values read from a file into a multi-array of size 3.
/*!
 * Function to parse a block of values read from a file into a multi-array of size 3.
 * The size of the independent variables is provided as input from this file.
 * \param independentVariableSize Size of independent variables (must be of size 3)
 * \param coefficientsBlock Raw number block as read from file.
 * \return Multi-array of size 3, according to input data.
 */
boost::multi_array< double, 3 > parseRawThreeDimensionalCoefficientsFromFile(
        const std::vector< int > independentVariableSize,
        const Eigen::MatrixXd& coefficientsBlock );

//! Function to parse a block of values read from a file into a multi-array of size 4.
/*!
 * Function to parse a block of values read from a file into a multi-array of size 4.
 * The size of the independent variables is provided as input from this file.
 * \param independentVariableSize Size of independent variables (must be of size 4)
 * \param coefficientsBlock Raw number block as read from file.
 * \return Multi-array of size 4, according to input data.
 */
boost::multi_array< double, 4 > parseRawFourDimensionalCoefficientsFromFile(
        const std::vector< int > independentVariableSize,
        const Eigen::MatrixXd& coefficientsBlock );

//! Function to read a coefficient file (data on a structured grid as a function of N independent variables)
/*!
 * Function to read a coefficient file (data on a structured grid as a function of N independent variables). This function
 * reads the file as raw data, converting the data into a multi-array of the correct size can be done using the
 * MultiArrayFileReader class if needed. The file format is defined in the Tudat wiki.
 * \param fileName Name of file containing coefficients,
 * \param independentVariables Independent variables of coefficients, as read from file (returned by reference)
 * \param coefficientBlock Block of coefficients, as read from file (returned by reference)
 */
void readCoefficientsFile(
        const std::string fileName,
        std::vector< std::vector< double > >& independentVariables,
        Eigen::MatrixXd& coefficientBlock );

//! Function to retrieve the number of independent variables that the coefficients in a file are given for.
/*!
 *  Function to retrieve the number of independent variables that the coefficients in a file are given for.
 *  \param fileName Name of the file containing the coefficients.
 *  \return Number of independent variables that the coefficients in a file are given for.
 */
int getNumberOfIndependentVariablesInCoefficientFile( const std::string& fileName );

//! Interface class for reading coefficients as a function of N independent variables from a file.
/*!
 *  Interface class for reading coefficients as a function of N independent variables from a file. This class is used instead
 *  of a single templated free function to  allow multi-arrays of different sizes to be created using the same interface.
 *  NOTE: The possibility of using a single  templated implementation for arbitrary multi-array size should be investigated
 *  in the future.
 */
template< unsigned int NumberOfDimensions >
class MultiArrayFileReader
{
public:
    
    //! Function to read only a multi-array from the file
    /*!
     *  Function to read only a multi-array from the file, not saving the values of the independent variables it is
     *  defined.
     *  \param fileName Name of the coefficient file
     *  \return Multi-array containing file contents.
     */
    static boost::multi_array< double, static_cast< size_t >( NumberOfDimensions ) > readMultiArray(
            const std::string fileName );
    
    //! Function to read a multi-array from the file, and the values of the independent variables at which it is defined
    /*!
     *  Function to read a multi-array from the file, and the values of the independent variables at which it is defined
     *  \param fileName Name of the coefficient file
     *  \return  Pair: first entry containing multi-array of double coefficients, second containing list of independent
     *  variables at which coefficients are defined.
     */
    static std::pair< boost::multi_array< double, 1 >, std::vector< std::vector< double > > >
    readMultiArrayAndIndependentVariables(
            const std::string fileName );

};

//! Interface class for reading coefficients as a function of 1 independent variables from a file.
template< >
class MultiArrayFileReader< 1 >
{
public:
    
    //! Function to read only a multi-array of size 1 from the file
    /*!
     *  Function to read only a multi-array of size 1 from the file, not saving the values of the independent variables it is
     *  defined.
     *  \param fileName Name of the coefficient file
     *  \return Multi-array containing file contents.
     */
    static boost::multi_array< double, 1 > readMultiArray(
            const std::string fileName )
    {
        return readMultiArrayAndIndependentVariables( fileName ).first;
    }
    
    //! Function to read a multi-array from the file, and the values of the independent variables at which it is defined
    /*!
     *  Function to read a multi-array of size 1 from the file, and the values of the independent variables at which it is
     *  defined
     *  \param fileName Name of the coefficient file
     *  \return  Pair: first entry containing multi-array of double coefficients, second containing list of independent
     *  variables at which coefficients are defined.
     */
    static std::pair< boost::multi_array< double, 1 >, std::vector< std::vector< double > > >
    readMultiArrayAndIndependentVariables(
            const std::string fileName )
    {
        std::vector< std::vector< double > > independentVariables;
        Eigen::MatrixXd coefficientBlock;
        
        // Read raw data from file
        readCoefficientsFile( fileName, independentVariables, coefficientBlock );
        
        if ( !( independentVariables.size( ) == 1 ) )
        {
            throw std::runtime_error( "Error when reading 3-D multi-array, wrong number of independent variables found"
                                      " (expected 1 independent variable, got " +
                                      std::to_string( independentVariables.size( ) ) + ")" );
        }
        
        // Set sizes of independent variables
        std::vector< int > independentVariableSize;
        independentVariableSize.push_back( independentVariables.at( 0 ).size( ) );
        
        // Parse data from file
        return std::make_pair( parseRawOneDimensionalCoefficientsFromFile(
                                   independentVariableSize, coefficientBlock ), independentVariables );
    }
};

//! Interface class for reading coefficients as a function of 2 independent variables from a file.
template< >
class MultiArrayFileReader< 2 >
{
public:
    
    //! Function to read only a multi-array of size 2 from the file
    /*!
     *  Function to read only a multi-array of size 2 from the file, not saving the values of the independent variables it is
     *  defined.
     *  \param fileName Name of the coefficient file
     *  \return Multi-array containing file contents.
     */
    static boost::multi_array< double, 2 > readMultiArray(
            const std::string fileName )
    {
        return readMultiArrayAndIndependentVariables( fileName ).first;
    }
    
    //! Function to read a multi-array from the file, and the values of the independent variables at which it is defined
    /*!
     *  Function to read a multi-array of size 2 from the file, and the values of the independent variables at which it is
     *  defined
     *  \param fileName Name of the coefficient file
     *  \return  Pair: first entry containing multi-array of double coefficients, second containing list of independent
     *  variables at which coefficients are defined.
     */
    static std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > >
    readMultiArrayAndIndependentVariables(
            const std::string fileName )
    {
        std::vector< std::vector< double > > independentVariables;
        Eigen::MatrixXd coefficientBlock;
        
        // Read raw data from file
        readCoefficientsFile( fileName, independentVariables, coefficientBlock );
        
        if ( !( independentVariables.size( ) == 2 ) )
        {
            throw std::runtime_error( "Error when reading 3-D multi-array, wrong number of independent variables found"
                                      " (expected 2 independent variables, got " +
                                      std::to_string( independentVariables.size( ) ) + ")" );
        }
        
        // Set sizes of independent variables
        std::vector< int > independentVariableSize;
        independentVariableSize.push_back( independentVariables.at( 0 ).size( ) );
        independentVariableSize.push_back( independentVariables.at( 1 ).size( ) );
        
        // Parse data from file
        return std::make_pair( parseRawTwoDimensionalCoefficientsFromFile(
                                   independentVariableSize, coefficientBlock ), independentVariables );
    }
};

//! Interface class for reading coefficients as a function of 3 independent variables from a file.
template< >
class MultiArrayFileReader< 3 >
{
public:
    
    //! Function to read only a multi-array of size 3 from the file
    /*!
     *  Function to read only a multi-array of size 3 from the file, not saving the values of the independent variables it is
     *  defined.
     *  \param fileName Name of the coefficient file
     *  \return Multi-array containing file contents.
     */
    static boost::multi_array< double, 3 > readMultiArray(
            const std::string fileName )
    {
        return readMultiArrayAndIndependentVariables( fileName ).first;
    }
    
    //! Function to read a multi-array from the file, and the values of the independent variables at which it is defined
    /*!
     *  Function to read a multi-array of size 3 from the file, and the values of the independent variables at which it is
     *  defined
     *  \param fileName Name of the coefficient file
     *  \return  Pair: first entry containing multi-array of double coefficients, second containing list of independent
     *  variables at which coefficients are defined.
     */
    static std::pair< boost::multi_array< double, 3 >, std::vector< std::vector< double > > >
    readMultiArrayAndIndependentVariables(
            const std::string fileName )
    {
        std::vector< std::vector< double > > independentVariables;
        Eigen::MatrixXd coefficientBlock;
        
        // Read raw data from file
        readCoefficientsFile( fileName, independentVariables, coefficientBlock );
        
        if ( !( independentVariables.size( ) == 3 ) )
        {
            throw std::runtime_error( "Error when reading 3-D multi-array, wrong number of independent variables found"
                                      " (expected 3 independent variables, got " +
                                      std::to_string( independentVariables.size( ) ) + ")" );
        }
        
        // Set sizes of independent variables
        std::vector< int > independentVariableSize;
        independentVariableSize.push_back( independentVariables.at( 0 ).size( ) );
        independentVariableSize.push_back( independentVariables.at( 1 ).size( ) );
        independentVariableSize.push_back( independentVariables.at( 2 ).size( ) );
        
        // Parse data from file
        return std::make_pair(
                    parseRawThreeDimensionalCoefficientsFromFile( independentVariableSize, coefficientBlock ),
                    independentVariables );
    }

};

//! Interface class for reading coefficients as a function of 3 independent variables from a file.
template< >
class MultiArrayFileReader< 4 >
{
public:

    //! Function to read only a multi-array of size 4 from the file
    /*!
     *  Function to read only a multi-array of size 4 from the file, not saving the values of the independent variables it is
     *  defined.
     *  \param fileName Name of the coefficient file
     *  \return Multi-array containing file contents.
     */
    static boost::multi_array< double, 4 > readMultiArray(
            const std::string fileName )
    {
        return readMultiArrayAndIndependentVariables( fileName ).first;
    }

    //! Function to read a multi-array from the file, and the values of the independent variables at which it is defined
    /*!
     *  Function to read a multi-array of size 3 from the file, and the values of the independent variables at which it is
     *  defined
     *  \param fileName Name of the coefficient file
     *  \return  Pair: first entry containing multi-array of double coefficients, second containing list of independent
     *  variables at which coefficients are defined.
     */
    static std::pair< boost::multi_array< double, 4 >, std::vector< std::vector< double > > >
    readMultiArrayAndIndependentVariables(
            const std::string fileName )
    {
        std::vector< std::vector< double > > independentVariables;
        Eigen::MatrixXd coefficientBlock;

        // Read raw data from file
        readCoefficientsFile( fileName, independentVariables, coefficientBlock );

        if ( !( independentVariables.size( ) == 4 ) )
        {
            throw std::runtime_error( "Error when reading 4-D multi-array, wrong number of independent variables found"
                                      " (expected 4 independent variables, got " +
                                      std::to_string( independentVariables.size( ) ) + ")" );
        }

        // Set sizes of independent variables
        std::vector< int > independentVariableSize;
        independentVariableSize.push_back( independentVariables.at( 0 ).size( ) );
        independentVariableSize.push_back( independentVariables.at( 1 ).size( ) );
        independentVariableSize.push_back( independentVariables.at( 2 ).size( ) );
        independentVariableSize.push_back( independentVariables.at( 3 ).size( ) );

        // Parse data from file
        return std::make_pair(
                    parseRawFourDimensionalCoefficientsFromFile( independentVariableSize, coefficientBlock ),
                    independentVariables );
    }

};

} // namespace input_output

} // namespace tudat

#endif // TUDAT_MULTIDIMENSIONALARRAYREADER_H

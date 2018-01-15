/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      The function printStandardScientificNotation() has been implemented to cope with
 *      cross-platform incompatibilities in the printed output of floating-point numbers in
 *      scientific notation.
 *
 */

#ifndef TUDAT_BASIC_INPUT_OUTPUT_H
#define TUDAT_BASIC_INPUT_OUTPUT_H

#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <Eigen/Core>

#include <boost/filesystem.hpp>

#include "Tudat/InputOutput/streamFilters.h"

namespace tudat
{
namespace input_output
{

//! Get root-path for Tudat library.
/*!
 * Returns root-path corresponding with root-directory of Tudat library as a string with
 * trailing slash included.
 * \return Tudat root-path.
 */
static inline std::string getTudatRootPath( )
{
#ifdef TUDAT_CUSTOM_ROOT_PATH
    return std::string( TUDAT_CUSTOM_ROOT_PATH );
#else
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path in the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    return filePath_.substr( 0, filePath_.length( ) -
                             std::string( "InputOutput/basicInputOutput.h" ).length( ) );
#endif
}

static inline std::string getEarthOrientationDataFilesPath( )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path in the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    return ( filePath_.substr( 0, filePath_.length( ) -
                               std::string( "InputOutput/basicInputOutput.h" ).length( ) ) +
             "Astrodynamics/EarthOrientation/Data/" );
}

//! Get Spice kernel path.
/*!
 * Returns path in which Spice kernels are located.
 * \return Path containing Spice kernels.
 */
static inline std::string getSpiceKernelPath( )
{
#ifdef SPICE_KERNEL_CUSTOM_FOLDER
    return std::string( SPICE_KERNEL_CUSTOM_FOLDER );
#else
    return getTudatRootPath( ) + "External/SpiceInterface/Kernels/";
#endif
}

//! Get atmosphere tables path.
/*!
 * Returns path in which atmosphere tables are located.
 * \return Path containing atmosphere tables.
 */
static inline std::string getAtmosphereTablesPath( )
{
    return getTudatRootPath( ) + "External/AtmosphereTables/";
}

//! Get gravity models path.
/*!
 * Returns path in which gravity models are located.
 * \return Path containing gravity models.
 */
static inline std::string getGravityModelsPath( )
{
    return getTudatRootPath( ) + "External/GravityModels/";
}

//! Get space weather data files path.
/*!
 * Returns path in which space weather data files are located.
 * \return Path containing space weather data files.
 */
static inline std::string getSpaceWeatherDataPath( )
{
    return getTudatRootPath( ) + "External/SpaceWeatherData/";
}

//! Print floating-point number in formatted scientific notation.
/*!
 * Prints floating-point number in formatted scientific notation. The user can specify the
 * precision that the base is represented by (default=maximum precision for doubles), and the
 * number of digits (width) in the exponent (default=2). This function can be used to ensure that
 * the string output generated in scientific notation is the consistently the same for all
 * platforms.
 * \param floatingPointNumber Floating-point number to format.
 * \param basePrecision Number of digits to represent the base of the floating-point number.
 * \param exponentWidth Number of digits to represent the exponent of the floating-point number.
 * \return Formatted string representation of floating-point number in scientific notation.
 */
std::string printInFormattedScientificNotation(
        const double floatingPointNumber,
        const int basePrecision = std::numeric_limits< double >::digits10,
        const int exponentWidth = 2 );

//! List all files in directory.
/*!
 * Lists all files in a given directory. There is a recursion option to allow
 * all files in subdirectories to be listed as well.
 * \param directory Absolute directory path.
 * \param isRecurseIntoSubdirectories Flag to set if algorithm should recurse through
 *          subdirectories. Set to false by default.
 * \return Container of filenames in directory, stored as Boost path variables.
 */
std::vector< boost::filesystem::path > listAllFilesInDirectory(
        const boost::filesystem::path& directory, const bool isRecurseIntoSubdirectories = false );

//! Write a value to a stream.
/*!
 * Write a value to a stream, left-aligned at a specified precision. Value is preceded by
 * delimiter followed by a space, and followed by an end-of-line character. Helper function for
 * writeDataMapToTextFile().
 * \tparam OutputStream Type of the stream to write to.
 * \tparam ValueType Type of the value to write, should support the << operator.
 * \param stream Stream to write to.
 * \param value Value to write to stream.
 * \param precision Precision to write the value with.
 * \param delimiter Delimiter to precede the value.
 */
template < typename OutputStream, typename ValueType >
void writeValueToStream( OutputStream& stream, const ValueType& value,
                         const int precision, const std::string& delimiter )
{
    stream << delimiter << " "
           << std::setprecision( precision ) << std::left
           << std::setw( precision + 1 ) << value << std::endl;
}

//! Write an Eigen type to a stream.
/*!
 * Write an Eigen type to a stream, row-by-row and left-aligned at a specified precision. Each
 * value is preceded by delimiter followed by a space, and followed by an end-of-line character.
 * Helper function for writeDataMapToTextFile().
 * \tparam OutputStream Type of the stream to write to.
 * \tparam ScalarType Scalar type for Eigen::Matrix to write to stream.
 * \tparam NumberOfRows Number of rows in Eigen::Matrix to write to stream.
 * \tparam NumberOfColumns Number of columns in Eigen::Matrix to write to stream.
 * \tparam Options options for Eigen::Matrix to write to stream.
 * \tparam MaxiumumRows Maximum number of rows for Eigen::Matrix to write to stream.
 * \tparam MinimumRows Minimum number of rows for Eigen::Matrix to write to stream.
 * \param stream Stream to write to.
 * \param value Value to write to stream.
 * \param precision Precision to write the value with.
 * \param delimiter Delimiter to precede the value.
 * \param endLineAfterRow Boolean to denote whether a new line is to be started after each row of the matrix
 */
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

//! Write data map to text file.
/*!
 * Writes data stored in a map to text file.
 * \tparam InputIterator Iterator of pairs.
 * \param iteratorDataMap Iterator to the first item of the sequence to write.
 * \param last Iterator to the item past the the last item of the sequence to write.
 * \param outputFilename Output filename.
 * \param outputDirectory Output directory. This can be passed as a string as well. It will be
 *          created if it does not exist.
 * \param fileHeader Text to be placed at the head of the output file. N.B: This string MUST end in
 *          a newline/return character, or the first line of data will not be printed on a new
 *          line.
 * \param precisionOfKeyType Number of significant digits of KeyType-data to output.
 * \param precisionOfValueType Number of significant digits of ValueType-data to output.
 * \param delimiter Delimiter character, to delimit data entries in file.
 */
template< typename InputIterator >
void writeDataMapToTextFile(
        InputIterator iteratorDataMap, InputIterator last,
        const std::string& outputFilename,
        const boost::filesystem::path& outputDirectory, const std::string& fileHeader,
        const int precisionOfKeyType, const int precisionOfValueType,
        const std::string& delimiter )
{
    // Check if output directory exists; create it if it doesn't.
    if ( !boost::filesystem::exists( outputDirectory ) )
    {
        boost::filesystem::create_directories( outputDirectory );
    }

    // Open output file.
    std::string outputDirectoryAndFilename = outputDirectory.string( ) + "/" + outputFilename;
    std::ofstream outputFile_( outputDirectoryAndFilename.c_str( ) );

    // Write file header to file.
    outputFile_ << fileHeader;

    // Loop over map of propagation history.
    for ( ; iteratorDataMap != last; iteratorDataMap++ )
    {
        // Print map data to output file.
        outputFile_ << std::setprecision( precisionOfKeyType )
                    << std::left << std::setw( precisionOfKeyType + 1 )
                    << iteratorDataMap->first;
        writeValueToStream( outputFile_, iteratorDataMap->second, precisionOfValueType,
                            delimiter );
    }

    // Close output file.
    outputFile_.close( );
}

//! Write data map to text file.
/*!
 * Writes data stored in a map to text file.
 * \tparam KeyType Data type for map key.
 * \tparam ValueType Data type for map value.
 * \param dataMap Map with data.
 * \param outputFilename Output filename.
 * \param outputDirectory Output directory. This can be passed as a string as well. It will be
 *          created if it does not exist.
 * \param fileHeader Text to be placed at the head of the output file. N.B: This string MUST end in
 *          a newline/return character, or the first line of data will not be printed on a new
 *          line.
 * \param precisionOfKeyType Number of significant digits of KeyType-data to output.
 * \param precisionOfValueType Number of significant digits of ValueType-data to output.
 * \param delimiter Delimiter character, to delimit data entries in file.
 */
template< typename KeyType, typename ValueType >
void writeDataMapToTextFile(
        const std::map< KeyType, ValueType >& dataMap, const std::string& outputFilename,
        const boost::filesystem::path& outputDirectory, const std::string& fileHeader = "",
        const int precisionOfKeyType = 16, const int precisionOfValueType = 16,
        const std::string& delimiter = "\t" )
{
    writeDataMapToTextFile( dataMap.begin( ), dataMap.end( ),
                            outputFilename, outputDirectory, fileHeader,
                            precisionOfKeyType, precisionOfValueType, delimiter );
}

//! Write data map to text file.
/*!
 * Writes data stored in a map to text file, using default KeyType-precision and
 * ValueType-precision (digits10 from "limits" standard library), output directory
 * (Tudat root-path), and delimiter (space).
 * Do not use this function for Eigen types. Instead, use the function that accepts a map as input
 * parameter below, or the function with the full argument list. Otherwise, your input will be
 * truncated.
 * \tparam InputIterator Iterator of pairs.
 * \param first Iterator to the first item of the sequence to write.
 * \param last Iterator to the item past the the last item of the sequence to write.
 * \param outputFilename Output filename.
 */
template< typename InputIterator >
void writeDataMapToTextFile( InputIterator first, InputIterator last,
                             const std::string& outputFilename )
{
    writeDataMapToTextFile( first, last,
                            outputFilename,
                            getTudatRootPath( ),
                            "",
                            std::numeric_limits< typename
                            InputIterator::value_type::first_type >::digits10,
                            std::numeric_limits< typename
                            InputIterator::value_type::second_type >::digits10,
                            " " );
}

//! Write data map to text file.
/*!
 * Writes data stored in a map to text file, using default KeyType-precision and
 * ValueType-precision (digits10 from "limits" standard library), output directory
 * (Tudat root-path), and delimiter (space).
 * \tparam KeyType Data type for map key.
 * \tparam ValueType Data type for map value.
 * \param dataMap Map with data.
 * \param outputFilename Output filename.
 */
template< typename KeyType, typename ValueType >
void writeDataMapToTextFile( const std::map< KeyType, ValueType >& dataMap,
                             const std::string& outputFilename )
{
    return writeDataMapToTextFile( dataMap.begin( ), dataMap.end( ),
                                   outputFilename,
                                   getTudatRootPath( ),
                                   "",
                                   std::numeric_limits< KeyType >::digits10,
                                   std::numeric_limits< ValueType >::digits10,
                                   " " );
}

//! Write data map to text file.
/*!
 * Writes Eigen data stored in a map to text file, using default KeyType-precision and
 * ValueType-precision (digits10 from "limits" standard library), output directory
 * (Tudat root-path), and delimiter (space).
 * \tparam KeyType Data type for map key.
 * \tparam ValueType Data type for Eigen::Matrix, used as map value.
 * \tparam NumberOfRows Number of rows in Eigen::Matrix, used as map value.
 * \tparam NumberOfColumns Number of columns in Eigen::Matrix, used as map value.
 * \tparam Options options for Eigen::Matrix, used as map value.
 * \tparam MaxiumumRows Maximum number of rows for Eigen::Matrix, used as map value.
 * \tparam MinimumRows Minimum number of rows for Eigen::Matrix, used as map value.
 * \param dataMap Map with data.
 * \param outputFilename Output filename.
 */
template< typename KeyType, typename ScalarType,
          int NumberOfRows, int NumberOfColumns, int Options, int MaximumRows, int MaximumCols >
void writeDataMapToTextFile( const std::map< KeyType, Eigen::Matrix< ScalarType,
                             NumberOfRows, NumberOfColumns, Options,
                             MaximumRows, MaximumCols > >& dataMap,
                             const std::string& outputFilename )
{
    return writeDataMapToTextFile( dataMap.begin( ), dataMap.end( ),
                                   outputFilename,
                                   getTudatRootPath( ),
                                   "",
                                   std::numeric_limits< KeyType >::digits10,
                                   std::numeric_limits< ScalarType >::digits10,
                                   " " );
}

//! Write data map to text file.
/*!
 * Writes Eigen data stored in a map to text file, using default KeyType-precision and
 * ValueType-precision (digits10 from "limits" standard library), output directory
 * (Tudat root-path), and delimiter (space).
 * \tparam KeyType Data type for map key.
 * \tparam ValueType Data type for Eigen::Matrix, used as map value.
 * \param dataMap Map with data.
 * \param outputPath Path of the output file.
 * \param fileHeader Text to be placed at the head of the output file. N.B: This string MUST end in
 *          a newline/return character, or the first line of data will not be printed on a new
 *          line.
 * \param precision Number of significant digits of KeyType-data and ValueType-data to output.
 */
template< typename KeyType, typename ValueType >
void writeDataMapToTextFile(
        const std::map< KeyType, ValueType >& dataMap,
        const boost::filesystem::path& outputPath,
        const std::string& fileHeader,
        const int precision )
{
    writeDataMapToTextFile( dataMap.begin( ), dataMap.end( ),
                            outputPath.filename( ).string( ), outputPath.parent_path( ),
                            fileHeader, precision, precision, " " );
}

//! Write Eigen matrix to text file.
/*!
 * Write Eigen matrix to text file.
 * \param matrixToWrite Matrix that is to be written to a file
 * \param outputFilename Output filename.
 * \param precisionOfMatrixEntries Number of significant digits of Matrix output
 * \param outputDirectory Output directory. It can be passed as a string as well. The directory be created if it does not exist.
 * \param delimiter Delimiter character, to delimit data entries in file.
 * \param header Header to be added in the first line. Header char (e.g. %, #), if any, must be included in this string.
 * Must end in line break to make data matrix start in the next line.
 */
template< typename ScalarType, int NumberOfRows, int NumberOfColumns >
void writeMatrixToFile( Eigen::Matrix< ScalarType, NumberOfRows, NumberOfColumns > matrixToWrite,
                        const std::string& outputFilename,
                        const int precisionOfMatrixEntries = 16,
                        const boost::filesystem::path& outputDirectory = getTudatRootPath( ),
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

//! Typedef for double-KeyType, double-ValueType map.
/*!
 * Typedef for double-KeyType, double-ValueType map.
 */
typedef std::map< double, double > DoubleKeyTypeDoubleValueTypeMap;

//! Typedef for double-KeyType, VectorXd-ValueType map.
/*!
 * Typedef for double-KeyType, VectorXd-ValueType map.
 */
typedef std::map< double, Eigen::VectorXd > DoubleKeyTypeVectorXdValueTypeMap;

//! Typedef for double-KeyType, Vector3d-ValueType map.
/*!
 * Typedef for double-KeyType, Vector3d-ValueType map.
 */
typedef std::map< double, Eigen::Vector3d > DoubleKeyTypeVector3dValueTypeMap;

//! Typedef for double-KeyType, MatrixXd-ValueType map.
/*!
 * Typedef for double-KeyType, MatrixXd-ValueType map.
 */
typedef std::map< double, Eigen::MatrixXd > DoubleKeyTypeMatrixXdValueTypeMap;

//! Typedef for double-KeyType, Matrix3d-ValueType map.
/*!
 * Typedef for double-KeyType, Matrix3d-ValueType map.
 */
typedef std::map< double, Eigen::Matrix3d > DoubleKeyTypeMatrix3dValueTypeMap;

//! Typedef for int-KeyType, double-ValueType map.
/*!
 * Typedef for int-KeyType, double-ValueType map.
 */
typedef std::map< int, double > IntKeyTypeDoubleValueTypeMap;

} // namespace input_output
} // namespace tudat

#endif // TUDAT_BASIC_INPUT_OUTPUT_H

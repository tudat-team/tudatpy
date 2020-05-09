/*    Copyright (c) 2010-2020, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TWO_LINE_ELEMENTS_FILE_READER_H
#define TUDAT_TWO_LINE_ELEMENTS_FILE_READER_H

#include <fstream>
#include <map>
#include <string>
#include <vector>

#include <memory>

namespace tudat
{
namespace input_output
{

//! Two-Line Elements (TLE) file reader class
/*!
 * Definition of TLE file reader class. Reads and parses a TLE data file (in the same formatting as
 * files from the CelesTrak website) into a map with the name of the satellite as the key and the
 * concatenated two lines with the elements as the value.
 *
 * Usage: create an instantiation of the class using the constructor (provide file name and absolute
 * path); then call the openFile( ) function to open the file. To process the file and store its content
 * in the map, call readAndStoreData( ). Finally, the map can then be retrieved by calling
 */
class TwoLineElementsFileReader
{
public:

    //! Constructor.
    /*!
     * Constructor from file name and absolute directory path.
     * @param fileName Name of the file including extension.
     * @param absoluteDirectoryPath String containing the path of the file to be read, including trailing (back)slash.
     * @param numberOfHeaderLines Number of header lines in the file (default: 0).
     * @param lineCommentChar Character that indicates that a line is a comment and should be skipped.
     */
    TwoLineElementsFileReader( std::string& fileName, std::string& absoluteDirectoryPath,
    		unsigned int numberOfHeaderLines = 0, std::string lineCommentChar = "#" )
        : numberOfHeaderLines_( numberOfHeaderLines ),
          dataFile_( ),
          fileName_( fileName ),
          stringOfData_( "" ),
          absoluteFilePath_( "" ),
          absoluteDirectoryPath_( absoluteDirectoryPath ),
          relativeDirectoryPath_( "" ),
          lineCommentChar_( lineCommentChar ),
          numberOfObjects_( 0 )
    { }

    //! Default destructor
    virtual ~TwoLineElementsFileReader( )
    {
    	if ( dataFile_.is_open( ) )
		{
    		dataFile_.close( );
		}
    }

    //! Set absolute directory path.
    /*!
     * Sets absolute path to directory containing data file. If this is set, the
     * relative path will be cleared.
     * \param absoluteDirectoryPath Absolute path to directory containing data file.
     */
    void setAbsoluteDirectoryPath( std::string& absoluteDirectoryPath )
    {
        relativeDirectoryPath_ = "";
        absoluteDirectoryPath_ = absoluteDirectoryPath;
    }

    //! Set relative directory path.
    /*!
     * Sets relative path to directory containing data file, with respect to
     * root directory of Tudat library. If this is set, the absolute path will be cleared.
     * \param relativeDirectoryPath Relative directory path.
     */
    void setRelativeDirectoryPath( std::string relativeDirectoryPath )
    {
        absoluteDirectoryPath_ = "";
        relativeDirectoryPath_ = relativeDirectoryPath;
    }

    //! Set file name.
    /*!
     * Sets file name of data file.
     * \param fileName File name.
     */
    void setFileName( std::string fileName )
    {
    	fileName_ = fileName;
    }

    //! Open data file.
    /*!
     * Opens data file.
     */
    void openFile( );

    //! Set number of lines for file header.
    /*!
     * Sets number of lines for file header. This function defines the number of
     * lines of the file header, starting from the beginning of the file.
     * This function cannot be used in combination with the skipLines( ),
     * skipLinesStartingWithCharacter( ), and skipLinesWithKeyword( ) functions.
     * \param numberOfHeaderLines Number of lines for file header.
     */
    void setNumberOfHeaderLines( unsigned int numberOfHeaderLines )
    {
        numberOfHeaderLines_ = numberOfHeaderLines;
    }

    //! Read and store data from data file.
    void readAndStoreData( );

    //! Get number of objects.
    /*!
     * \return Number of objects in TLE data file.
     */
    unsigned int getNumberOfObjects( ) { return numberOfObjects_; }

    std::map< std::string, std::string > getTleMap( ) { return tlePairs_; };

protected:


    //! Number of header lines.
    unsigned int numberOfHeaderLines_;

    //! Data file stream.
    std::ifstream dataFile_;

    //! File name.
    std::string fileName_;

    //! String of data.
    std::string stringOfData_;

    //! Absolute path to data file.
    std::string absoluteFilePath_;

    //! Absolute path to directory containing data file.
    std::string absoluteDirectoryPath_;

    //! Relative path to directory containing data file.
    std::string relativeDirectoryPath_;

    //! Starting character used by the skipLinesStartingWithCharacter( ) function.
    std::string lineCommentChar_;

    std::map< std::string, std::string > tlePairs_;

    //! Container of header data from file.
    /*!
     * Map container of string header data from data file, obtained by reading each
     * line of header data from file using the getline( ) function. The key is the
     * line number from the file and the value is header line data.
     */
    std::vector< std::string > headerLines_;

private:

    //! Number of object in input file, linecounter divided by 3.
    unsigned int numberOfObjects_;

};

} // namespace input_output
} // namespace tudat

#endif // TUDAT_TWO_LINE_ELEMENTS_FILE_READER_H

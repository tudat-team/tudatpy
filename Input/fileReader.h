/*! \file fileReader.h
 *    This header file contains the definition of a base class for all file
 *    readers included in Tudat.
 *
 *    Path              : /Input/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J. Leloux
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.Leloux@student.tudelft.nl
 *
 *    Date created      : 24 February, 2011
 *    Last modified     : 27 June, 2011
 *
 *    References
 *
 *    Notes
 *      This class implements two ways to read input files. One makes use of
 *      the readAndStoreData( ) function with the
 *      skipLinesStartingWithCharacter( ) function, which results in an entire
 *      file being read and stored, except for the lines starting with the
 *      given starting character. The second way to read input files is by
 *      using the readAndStoreData( numberOfLines ) and
 *      skipLines( numberOfLines ) functions. This method does not skip lines
 *      starting with specific characters.
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110223    K. Kumar          First creation of code.
 *      110224    K. Kumar          Changed vector container to map container.
 *      110224    J. Leloux         Checked code and fixed typo.
 *      110627    K. Kumar          Moved skipLinesWithKeyword() from
 *                                  TextFileReader.
 */

#ifndef FILEREADER_H
#define FILEREADER_H

// Include statements.
#include <iostream>
#include <fstream>
#include <map>
#include "Basics/basicFunctions.h"

//! File reader class.
/*!
 * Definition of file reader class.
 */
class FileReader
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    FileReader( ) : lineCounter_( 1 ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~FileReader( ) { }

    //! Set relative path.
    /*!
     * Sets relative path to directory containing data file, with respect to
     * root directory of Tudat library.
     * \param relativePath Relative path.
     */
    void setRelativePath( std::string relativePath ) { relativePath_ = relativePath; }

    //! Set file name.
    /*!
     * Sets file name of data file.
     * \param fileName File name.
     */
    void setFileName( std::string fileName ) { fileName_ = fileName; }

    //! Open data file.
    /*!
     * Opens data file.
     */
    void openFile( );

    //! Skip lines.
    /*!
     * Skips a given number of lines of data file. This function cannot be used
     * in combination with the skipLinesStartingWithCharacter( ) function.
     * \param numberOfLines Number of lines.
     */
    void skipLines( unsigned int numberOfLines );

    //! Skip all lines starting with a given character.
    /*!
     * Skips all lines starting with a given character. This function cannot be
     * used in combination with the skipLines( ) function.
     * \param startingCharacter Starting character.
     */
    void skipLinesStartingWithCharacter( const std::string& startingCharacter )
    { startingCharacter_ = startingCharacter; }

    //! Skip all lines containing a given keyword.
    /*!
     * Skips all lines containing a given keyword. This function cannot be
     * used in combination with the skipLines( ) function.
     * \param skipKeyword Keyword to skip line.
     */
    void skipLinesWithKeyword( const std::string& skipKeyword ) { skipKeyword_ = skipKeyword; }

    //! Read and store data.
    /*!
     * Reads and stores data from data file. When used in combination with the
     * skipLinesStartingWithCharacter( ) function, this function will read all
     * lines of the input file, except those starting with the given starting
     * character.
     */
    virtual void readAndStoreData( ) = 0;

    //! Read and store data.
    /*!
     * Reads and stores given number of lines of data from data file. This
     * function cannot be used in combination with the
     * skipLinesStartingWithCharacter( ) function.
     * \param numberOfLines Number of lines.
     */
    virtual void readAndStoreData( unsigned int numberOfLines ) = 0;

    //! Close data file.
    /*!
     * Closes data file.
     */
    void closeFile( ) { dataFile_.close( ); }

    //! Get vector container of data from file.
    /*!
     * Returns map container of string data from data file.
     * \return Pointer to map container of data from file.
     */
    std::map< unsigned int, std::string >& getContainerOfData( )
    { return containerOfDataFromFile_; }

protected:

    //! Line counter.
    /*!
     * Line counter.
     */
    unsigned int lineCounter_;

    //! Data file stream.
    /*!
     * Data file stream.
     */
    std::ifstream dataFile_;

    //! File name.
    /*!
     * File name for data file.
     */
    std::string fileName_;

    //! String of data.
    /*!
     * String of data.
     */
    std::string stringOfData_;

    //! Absolute path to data file.
    /*!
     * Absolute path to data file.
     */
    std::string absolutePath_;

    //! Relative path to data file.
    /*!
     * Relative path to data file with respect to root directory of Tudat.
     */
    std::string relativePath_;

    //! Starting character.
    /*!
     * Starting character.
     */
    std::string startingCharacter_;

    //! Skip keyword.
    /*!
     * Keyword used to skip a line.
     */
    std::string skipKeyword_;

    //! Map container of data from file.
    /*!
     * Map container of string data from data file, obtained by reading each
     * line of data from file using the getline( ) function. The key is the
     * line number from the file and the value is line data.
     */
    std::map< unsigned int, std::string > containerOfDataFromFile_;

    //! Iterator for map container of data from file.
    /*!
     * Iterator for map container of string data from data file.
     */
    std::map< unsigned int, std::string >::iterator iteratorContainerOfDataFromFile_;

private:
};

#endif // FILEREADER_H

// End of file.

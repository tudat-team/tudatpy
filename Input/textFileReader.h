/*! \file textFileReader.h
 *    This header file contains the definition of a text file reader class.
 *
 *    Path              : /Input/
 *    Version           : 5
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
 *    Date created      : 23 February, 2011
 *    Last modified     : 27 June, 2011
 *
 *    References
 *
 *    Notes
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
 *      110315    J. Leloux         Checked code.
 *      110316    K. Kumar          Added ostream >> operator overload.
 *      110607    F.M. Engelen      Added skipKeyWord Feature.
 *      110627    K. Kumar          Moved skipLinesWithKeyword() to FileReader.
 */

#ifndef TEXTFILEREADER_H
#define TEXTFILEREADER_H

// Include statements.
#include <iostream>
#include "Input/fileReader.h"

//! Text file reader class.
/*!
 * Definition of text file reader class.
 */
class TextFileReader : public FileReader
{
public:

    //! Read and store data.
    /*!
     * Reads and stores data from data file.
     */
    void readAndStoreData( );

    //! Read and store data.
    /*!
     * Reads and stores given number of lines of data from data file.
     * \param numberOfLines Number of lines.
     */
    void readAndStoreData( unsigned int numberOfLines );

    //! Strip End-Of-Line characters.
    /*!
     * Strips End-Of-Line (EOL) characters from string data stored in
     * containerOfDataFromFile_. The EOL characters removed are "\r" and "\n".
     * This function should only be used once data input data from the assigned
     * data file has been read and stored.
     */
    void stripEndOfLineCharacters( );

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param pointerToTextFileReader Pointer to text file reader.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream,
                                     TextFileReader& pointerToTextFileReader );

protected:

private:
};

#endif // TEXTFILEREADER_H

// End of file.

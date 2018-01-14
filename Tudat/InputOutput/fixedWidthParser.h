/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_FIXED_WIDTH_PARSER_H
#define TUDAT_FIXED_WIDTH_PARSER_H

#include "Tudat/InputOutput/textParser.h"

namespace tudat
{
namespace input_output
{

//! Fixed-width parser.
/*!
 * This class can be used to parse a raw line of data into a series of fields, based on a fixed
 * field width.
 * This class inherits from the TextParser abstract base class.
 * \sa TextParser
 */
class FixedWidthParser : public TextParser
{
public:

    //! Create a parser that parses based on specified field widths and field type list.
    /*!
     * \param numberOfFields Number of fields to parse.
     * Arguments are the field types:
     * e.g.
     *     FixedWidthParser(4, fieldtypes::ID, fieldtypes::X, fieldtypes::Y, fieldtypes::Z);
     */
    FixedWidthParser( int numberOfFields, ... );

    //! Set trim: Trim whitespace off fields (default=true).
    void setTrim( bool trim ) { doTrim = trim; }

    //! Get trim setting: Trim whitespace off fields (default=true).
    bool getTrim( ) { return doTrim; }

    //! Set unit transformation map.
    void setUnitTransformationMap (
            std::map< FieldType, boost::shared_ptr< FieldTransform > > unitTransformationMap )
    {
        unitTransformationMap_ = unitTransformationMap;
    }

protected:

    //! Parses one line of text.
    /*!
     * Parses one line of text by dividing the line in fields using specified field widths.
     * Any white space at the start or end of a field are trimmed.
     * \param line Raw line of data.
     */
    void parseLine( std::string& line );

private:

    //! Number of fields that is parsed.
    unsigned int numberOfFields_;

    //! Vector containing FieldTypes of parsed data.
    std::vector< FieldType > typeList;

    //! Vector containing the widths of each field.
    std::vector< int > sizeList;

    //! Boolean: true = trim whitespace off fields; false = no trim.
    bool doTrim;

    //! Unit transformation map.
    std::map< FieldType, FieldTransformPointer > unitTransformationMap_;
};

//! Typedef for shared-pointer to FixedWidthParser object.
typedef boost::shared_ptr< FixedWidthParser > FixedWidthParserPointer;

} // namespace input_output
} // namespace tudat

#endif // TUDAT_FIXED_WIDTH_PARSER_H

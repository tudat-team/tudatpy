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

#ifndef TUDAT_SEPARATED_PARSER_H
#define TUDAT_SEPARATED_PARSER_H

#include <cstdarg>
#include <string>

#include <memory>

#include "Tudat/InputOutput/textParser.h"

namespace tudat
{
namespace input_output
{

//! Separated parser class.
/*!
 * The separated parsed class can be used to parse a raw line of data into a series of values,
 * separated by a specified separator string, e.g., a space, a comma, a tab etc.
 * This class inherits from the TextParser abstract base class.
 * \sa TextParser
 */
class SeparatedParser : public TextParser
{
public:

    //! Create a parser that parses based on a specified separator and field type list.
    /*!
     * \param separator String representation of symbol or text which is used as separator.
     * \param numberOfFields Number of fields to parse.
     * Arguments are the field types:
     * e.g.
     *     SeparatedParser(",", 4, fieldtypes::ID, fieldtypes::X, fieldtypes::Y, fieldtypes::Z);
     */
    SeparatedParser( std::string separator, int numberOfFields, ... );

    //! Set trim: Trim whitespace off fields (default=true).
    void setTrim( bool trim ) { doTrim = trim; }

    //! Get trim setting: Trim whitespace off fields (default=true).
    bool getTrim( ) { return doTrim; }

    //! Set unit transformation map.
    /*!
     * \param unitTransformationMap Map providing field transforms for any or all of the
     *          FieldTypes.
     */
    void setUnitTransformationMap(
            std::map< FieldType, std::shared_ptr< FieldTransform > > unitTransformationMap )
    {
        unitTransformationMap_ = unitTransformationMap;
    }

protected:

    //! Parses one line of text.
    /*!
     * Parses one line of text by dividing the line in fields using the specified separator.
     * \param line Raw line of data.
     */
    void parseLine( std::string& line );

private:

    //! Number of fields that is parsed.
    unsigned int numberOfFields_;

    //! String containing the separator.
    std::string separator_;

    //! Vector containing FieldTypes of parsed data.
    std::vector< FieldType > typeList;

    //! Boolean: true = trim whitespace off fields; false = no trim.
    bool doTrim;

    //! Map containing unit transformation equations.
    /*!
     * Map containing equations which are used to transform raw data to SI unit data.
     */
    std::map< FieldType, FieldTransformPointer > unitTransformationMap_;
};

//! Typedef for shared-pointer to SeparatedParser object.
typedef std::shared_ptr< SeparatedParser > SeparatedParserPointer;

} // namespace input_output
} // namespace tudat

#endif // TUDAT_SEPARATED_PARSER_H

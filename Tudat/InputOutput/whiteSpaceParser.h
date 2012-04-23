/*!   Copyright (c) 2010-2011 Delft University of Technology.
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
 *      111103    S. Billemont      First creation of code.
 *      120217    D.J. Gondelach    Code check.
 */

#ifndef TUDAT_WHITESPACEPARSER_H
#define TUDAT_WHITESPACEPARSER_H

#include <string>
#include <stdarg.h>

#include <boost/shared_ptr.hpp>

#include "Tudat/InputOutput/textParser.h"

namespace tudat
{
namespace input_output
{

//! Class that parses lines based on white spaces.
/*!
 * A WhiteSpaceParser is a text parser that parses lines based on white spaces. The input string
 * is splitted into single data values (field values) which are associated with their field types.
 * There is no limit on the number of white spaces before, between and after data entries. However,
 * tabs will cause an error! This parser verifies that the number of field values found corresponds
 * to the number of specified field types. Field types should be specified in order of appearance
 * along the data string.
 *
 * NOTE: This WhiteSpaceParser works with the FieldValue/FieldType architecture.
 * For simpler file reading, use, for instance, matrixTextFileReader.
 */
class WhiteSpaceParser : public TextParser
{

public:

    //! Create a parser that parses based on white spaces given a field type list.
    /*!
     * \param numberOfFields Number fields to parse per line
     * Arguments passed as ... are the field types:
     * e.g.
     *      WhiteSpaceParser(4, fieldtypes::ID, fieldtypes::X, fieldtypes::Y, fieldtypes::Z);
     */
    WhiteSpaceParser( int numberOfFields, ... );

    //! Default destructor.
    virtual ~WhiteSpaceParser( ) { }

    //! Set unit transformation map.
    /*!
     * \param unitTransformationMap Map providing field transforms for any or all of the FieldTypes.
     */
    void setUnitTransformationMap(
            std::map< FieldType, boost::shared_ptr< FieldTransform > > unitTransformationMap )
    {
        unitTransformationMap_ = unitTransformationMap;
    }

protected:

    //! Parse line of data.
    /*!
     *  Parse a line of data by distinguish data fields using white spaces.
     * \param line String to parse
     */
    void parseLine( std::string& line );

private:

    //! Number of fields that is parsed.
    unsigned int numberOfFields_;

    //! Vector containing FieldTypes of parsed data.
    std::vector< FieldType > typeList;

    //! Map containing unit transformation equations.
    /*!
     *  Map containing equations which are used to transform raw data to SI unit data.
     */
    std::map< FieldType, boost::shared_ptr<FieldTransform > > unitTransformationMap_;

};

} // namespace input_output
} // namespace tudat
#endif

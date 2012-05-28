/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      111103    S. Billemont      Creation of code.
 *      120217    D.J. Gondelach    Code check.
 *
 *    References
 *
 */

#ifndef TUDAT_WHITE_SPACE_PARSER_H
#define TUDAT_WHITE_SPACE_PARSER_H

#include <cstdarg>
#include <map>
#include <string>
#include <vector>

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
    /*!
     * Default destructor.
     */
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
     * Map containing equations which are used to transform raw data to SI unit data.
     */
    std::map< FieldType, boost::shared_ptr<FieldTransform > > unitTransformationMap_;

};

} // namespace input_output
} // namespace tudat

#endif // TUDAT_WHITE_SPACE_PARSER_H

/*    Copyright (c) 2010-2013, Delft University of Technology
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
 *      111103    S. Billemont      File created.
 *      120706    D.J. Gondelach    Code check.
 *      130121    K. Kumar          Added shared-ptr typedef.
 *      131221    K. Kumar          Fixed missing Doxygen comments.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_SEPARATED_PARSER_H
#define TUDAT_SEPARATED_PARSER_H

#include <cstdarg>
#include <string>

#include <boost/shared_ptr.hpp>

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
            std::map< FieldType, boost::shared_ptr< FieldTransform > > unitTransformationMap )
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
typedef boost::shared_ptr< SeparatedParser > SeparatedParserPointer;

} // namespace input_output
} // namespace tudat

#endif // TUDAT_SEPARATED_PARSER_H

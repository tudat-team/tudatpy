/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      120326    D. Dirkx          Code checked, minor layout changes, implementation moved
 *                                  to cpp file and static identifier removed to prevent compile
 *                                  warning.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_PARSED_DATA_VECTOR_UTILITIES_H
#define TUDAT_PARSED_DATA_VECTOR_UTILITIES_H

#include <cstdarg>
#include <iostream>
#include <map>
#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/InputOutput/fieldType.h"
#include "Tudat/InputOutput/fieldValue.h"

namespace tudat
{
namespace input_output
{
namespace parsed_data_vector_utilities
{

//! Pointer to a field value (string).
typedef boost::shared_ptr< FieldValue >         FieldValuePtr;

//! Map containing field value pointers (mapped value), identified by their field type (key
//! value). Such maps contain the information from one parsed line.
typedef std::map< FieldType, FieldValuePtr >    ParsedDataLineMap;

//! Pointer to a parsed data line (see ParsedDataLineMap).
typedef boost::shared_ptr< ParsedDataLineMap >  ParsedDataLineMapPtr;

//! Vector of data map pointers (see ParsedDataLineMapPtr). Each entry of the vector points to
//! parsed data line.
typedef std::vector< ParsedDataLineMapPtr >     ParsedDataVector;

//! Pointer to the data map pointer vectors (see ParsedDataVector).
typedef boost::shared_ptr< ParsedDataVector >   ParsedDataVectorPtr;

//! Get the value of a given field from the map containing data.
/*!
 * Returns the transformed ( calls FieldValue::get() ) value of a field in a data map.
 *
 * Note: You must be sure that that the key exists in the data map or otherwise, segfaults will
 * occur.
 *
 * Example usage:
 *   ParsedDataLineMapPtr datamap;
 *   string planetName         = getField<string>(datamap, fieldtypes::general::name);
 *   int    planetID           = getField< int  >(datamap, fieldtypes::general::id);
 *   double planetEccentricity = getField<double>(datamap, fieldtypes::state::eccentricity);
 *
 * \param data   Check this datamap for all the requested FieldTypes.
 * \param field  Fieldtype to search for.
 * \return FieldValue The FieldValue string converted to the given type.
 */
template< typename T >
inline T getField( ParsedDataLineMapPtr data, FieldType field )
{
    boost::shared_ptr< FieldValue > fieldValue = data->find( field )->second;
    return fieldValue->get< T >( );
}

//! Get a shared-pointer to the value of a given field from the map containing data.
/*!
 * Returns the transformed ( calls FieldValue::getPointer() ) value of a field in a data map.
 *
 * Note: You must be sure that that the key exists in the data map or otherwise, segfaults will
 * occur.
 *
 * Example usage:
 *   ParsedDataLineMapPtr datamap;
 *   boost::shared_ptr<string> planetName
 *                                    = getField<string>(datamap, fieldtypes::general::name);
 *   boost::shared_ptr<int>    planetID
 *                                    = getField< int  >(datamap, fieldtypes::general::id);
 *   boost::shared_ptr<double> planetEccentricity
 *                                    = getField<double>(datamap, fieldtypes::state::eccentricity);
 *
 * \param data   Check this datamap for all the requested FieldTypes.
 * \param field  Fieldtype to search for.
 * \return FieldValue The FieldValue string converted to the given type.
 */
template< typename T >
inline boost::shared_ptr< T > getFieldPointer( ParsedDataLineMapPtr data, FieldType field )
{
    boost::shared_ptr< FieldValue > fieldValue = data->find( field )->second;
    return fieldValue->getPointer< T >( );
}

//! Filter the data vector for entries containing a given FieldType.
/*!
 * This allows for filtering of a parsed data vector, based on the requirement that each of the
 * datamaps MUST include each passed FieldTypes.
 *
 * Example usage:
 *   ParsedDataVectorPtr datavector;
 *   filterMapKey(datavector, 2, fieldtypes::general::name, fieldtypes::time::epoch);
 *
 * This return a new data vector that contains data maps only with the name
 * FieldType and the epoch FieldType.
 *
 * \param datavector Data vector to be checked for all the requested FieldTypes.
 * \param nrFields   Number of field that each entry must contain
 * \param ...        Fieldtypes that each entry must contain.
 * \return newdatavector A copy of the data vector with only the entries that contain the
 *                       passed FieldTypes.
 */
ParsedDataVectorPtr filterMapKey( ParsedDataVectorPtr datavector, int nrFields, ... );

//! Filter the data vector vector for entries containing a given FieldType and a matching
//! FieldValue regex.
/*!
 * This allows for filtering of a parsed data vector, based on the requirement that each of the
 * datamaps MUST include each passed FieldType and its value MUST have at least one match for
 * the passed regex for each respective field.
 *
 * Example usage:
 *   ParsedDataVectorPtr datavector;
 *   filterMapKeyValue(datavector, 2, fieldtypes::general::name, "Venus",
 *                                    fieldtypes::time::epoch, "2451545[\.0-9]*");
 *
 * This filter would create a new vector, where each data line is about Venus, and the epoch
 * for that dataline is on the date of 2000-1-1 (J2000). This means epoch can range from
 * 2000-1-1 00h to 2000-1-1 23h59m59.999...
 *
 * Note: the regular expression is in the form of char* and not boost::regex because that
 * yields the followong stdargs incompatibility:
 * error: cannot receive objects of non-trivially-copyable type 'boost::regex {aka struct
 * boost::basic_regex<char, boost::regex_traits<char> >}' through '...';
 *
 * \param datavector Data vector to be checked for all the requested FieldTypes.
 * \param nrFields   Number of field that each entry must contain.
 * \param ...        List of FieltType and Regex (char*) entries that must match to be addmited
 *                   in the filtered vector.
 * \return newdatavector    A copy of the data vector with only the entries that contain the
 *                          passed FieldTypes.
 */
ParsedDataVectorPtr filterMapKeyValue( ParsedDataVectorPtr datavector, int nrFields, ... );

//! Dump the content of a data map to an ostream.
/*!
 * This method allows users to quickly visualize the content of a datamap (single line)
 * by dumping its raw or transformed content to a given output stream (eg std::cout).
 *
 * Example usage:
 *   ParsedDataLineMapPtr datamap;
 *   dump(std::cout, datamap, false);
 *
 * This dumps untransformed content of datamap to std::cout
 *
 * \param stream            Stream to dump the content to.
 * \param data              Datamap content to dump.
 * \param showTransformed   Flag to show if it should dump raw values (false)
 *                          or transformed values (true).
 * \return stream A reference to the ostream.
 */
std::ostream& dump( std::ostream& stream, ParsedDataLineMapPtr data, bool showTransformed );

//! Dump the content of a data vector to an ostream (eg std::cout)
/*!
 * This method allows users to quickly visualize the content of a data vector by dumping
 * its raw or transformed content to a given output stream.
 *
 * Example usage:
 *   ParsedDataVectorPtr datavector;
 *   dump(std::cout, datavector, false);
 *
 * This dump untransformed content of datavector to std::cout.
 *
 * \param stream            Stream to dump the content to.
 * \param data              Datavector content to dump.
 * \param showTransformed   Flag to show if it should dump raw values (false)
 *                           or transformed values (true)
 * \return A reference to the ostream
 */
std::ostream& dump( std::ostream& stream, ParsedDataVectorPtr data, bool showTransformed );

} // namespace parsed_data_vector_utilities
} // namespace input_output
} // namespace tudat

#endif // TUDAT_PARSED_DATA_VECTOR_UTILS_H

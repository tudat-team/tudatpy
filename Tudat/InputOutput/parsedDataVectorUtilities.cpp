/* Copyright (c) 2010-2011 Delft University of Technology.
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
 *      120326    D. Dirkx          Code checked, minor layout changes, implementation moved
 *                                  to cpp file and static identifier removed to prevent compile
 *                                  warning.
 */

#include "Tudat/InputOutput/parsedDataVectorUtilities.h"

namespace tudat
{
namespace input_output
{
namespace parsed_data_vector_utilities
{

    //! Get the value of a given field from the data map containing (K=FieldType, V=fieldvalue).
    template< typename V >
     V getField( ParsedDataLineMapPtr data, FieldType field )
    {
        boost::shared_ptr<FieldValue> value = data->find(field)->second;
        boost::shared_ptr<std::string> str = value->get();
        return boost::lexical_cast<V>(*str);
    }

    //! Filter the data vector for entries containing a given FieldType.
     ParsedDataVectorPtr filterMapKey(ParsedDataVectorPtr datavector, int nrFields, ...)
    {
        // Create a fancy vector (list) of all the fields:
        va_list	argumentList;               // Define argument list variable.
        va_start(argumentList, nrFields);   // Initialize list; point to last defined argument.

        // Create a new datavector for the filtered data.
        ParsedDataVectorPtr newdatavector = ParsedDataVectorPtr(new ParsedDataVector());

        // Make a simple list to iterate over from all the FieldType arguments.
        std::vector<FieldType> checkForFieldTypes;
        for (int i=0; i < nrFields; i++)
        {
            checkForFieldTypes.push_back(va_arg(argumentList, FieldType));
        }

        // Clean up the system stack.
        va_end (argumentList);

        // Go over every dataline in the current datavector.
        for (ParsedDataVector::iterator currentDataLine = datavector->begin();
             currentDataLine != datavector->end();
             currentDataLine++)
        {
            // Flag to indicate that all the FieldTypes from checkForFieldTypes are present in this
            // line.
            bool found = true;

            // Loop over each FieldType to check if it exists.
            for(std::vector<FieldType>::iterator currentFieldCheck = checkForFieldTypes.begin();
                currentFieldCheck != checkForFieldTypes.end();
                currentFieldCheck++)
            {
                // Check if the FieldType is in the current data line.
                if ((*currentDataLine)->find(*currentFieldCheck) == (*currentDataLine)->end())
                {
                    // If not, mark that not all entries are present, and stop searching for others.
                    found = false;
                    break;
                }

                // If all the fields are present, add the current line to the new (filtered) vector.
                if (found)
                    newdatavector->push_back(*currentDataLine);
            }
        }

        // Return the filtered data vector
        return newdatavector;
    }

    //! Filter the data vector vector for entries containing a given FieldType and a matching
     ParsedDataVectorPtr filterMapKeyValue(ParsedDataVectorPtr datavector, int nrFields, ...)
    {
        // Create a fancy vector (list) of all the fields:
        va_list	argumentList;               // Define argument list variable.
        va_start(argumentList, nrFields);   // Initialize list; point to last defined argument.

        // Create a new data vector for the filtered data.
        ParsedDataVectorPtr newdatavector = ParsedDataVectorPtr(new ParsedDataVector());

        // Make a simple list to iterate over with the FieldType arguments and respective regex
        // expressions.
        std::map<FieldType, boost::regex> checkForFieldTypes;
        for (int i=0; i < nrFields; i++)
        {
            FieldType type = va_arg(argumentList, FieldType);
            boost::regex regex = boost::regex(va_arg(argumentList, char*));
            checkForFieldTypes.insert(std::pair<FieldType, boost::regex>(type, regex));
        }

        // Clean up the system stack.
        va_end (argumentList);

        // Go over every dataline in the current data vector.
        for (ParsedDataVector::iterator currentDataLine = datavector->begin();
             currentDataLine != datavector->end();
             currentDataLine++)
        {
            // Flag to indicate that all the FieldTypes from checkForFieldTypes are present in this
            // line.
            bool found = true;

            // Loop over each FieldType to check if it exists.
            for(std::map<FieldType, boost::regex>::iterator currentFieldCheck
                = checkForFieldTypes.begin();
                currentFieldCheck != checkForFieldTypes.end();
                currentFieldCheck++)
            {
                // Check if the FieldType is in the current dataline
                ParsedDataLineMap::iterator entry =
                        (*currentDataLine)->find(currentFieldCheck->first);
                if (entry == (*currentDataLine)->end())
                {
                    // If not, mark that not all entry pairs are present, and stop searching for
                    // others.
                    found = false;
                    break;
                }

                // Get the field value string.
                boost::shared_ptr<std::string> str = entry->second->get();

                // Check if the field value regex matches.
                if (!boost::regex_search(str->begin(), str->end(), currentFieldCheck->second))
                {
                    // If not, mark this, and stop searching for others.
                    found = false;
                    break;
                }
            }

            // If all the fields are present, add the current line to the new (filtered) vector.
            if (found)
                newdatavector->push_back(*currentDataLine);
        }

        // Return the filtered data vector.
        return newdatavector;
    }

    //! Dump the content of a data map to an ostream.
     std::ostream& dump(std::ostream& stream, ParsedDataLineMapPtr data, bool showTransformed)
    {
        // Initial character to separate elements.
        std::cout << "|";

        // Loop over data map.
        for(ParsedDataLineMap::iterator element = data->begin(); element != data->end(); element++)
        {
            // If transformed flag is true, dump transformed values.
            if (showTransformed)
                stream << element->second->get()->c_str();
            // If not, dump raw values.
            else
                stream << element->second->getRaw()->c_str();

            // Final character to separate elements.
            stream << "\t|";
        }

        stream << std::endl;
        return stream;
    }

    //! Dump the content of a data vector to an ostream (eg std::cout).
    std::ostream& dump(std::ostream& stream, ParsedDataVectorPtr data,
                             bool showTransformed)
    {
        // Loop over vector.
        for(std::size_t i=0; i < data->size(); i++)
        {
            // Get data line.
            ParsedDataLineMapPtr line = data->at(i);

            // Call dump function for single line (datamap).
            stream << dump(stream, line, showTransformed);
        }
        return stream;
    }

} // namespace parsed_data_vector_utilities
} // namespace input_output
} // namespace tudat




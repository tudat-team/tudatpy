/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_OPTIONS_H
#define TUDAT_JSONINTERFACE_OPTIONS_H

#include "valueAccess.h"
#include "valueConversions.h"

namespace tudat
{

namespace json_interface
{

// ExceptionResponseType

//! Map of `ExceptionResponseType` string representations.
static std::map< ExceptionResponseType, std::string > exceptionResponseTypes =
{
    { continueSilently, "continueSilently" },
    { printWarning, "printWarning" },
    { throwError, "throwError" }
};

//! Convert `ExceptionResponseType` to `json`.
inline void to_json( json& jsonObject, const ExceptionResponseType& exceptionResponseType )
{
    jsonObject = json_interface::stringFromEnum( exceptionResponseType, exceptionResponseTypes );
}

//! Convert `json` to `ExceptionResponseType`.
inline void from_json( const json& jsonObject, ExceptionResponseType& exceptionResponseType )
{
    exceptionResponseType = json_interface::enumFromString( jsonObject.get< std::string >( ), exceptionResponseTypes );
}


// ApplicationOptions

//! Class containing the application options for Tudat apps that use json_interface.
class ApplicationOptions
{
public:
    //! Constructor.
    ApplicationOptions( ) { }

    //! Response to a "default value used for missing key" event.
    ExceptionResponseType defaultValueUsedForMissingKey = continueSilently;

    //! Response to a "unused key" event.
    ExceptionResponseType unusedKey = printWarning;

    //! Response to a "unidimensional array inference" event.
    ExceptionResponseType unidimensionalArrayInference = continueSilently;

    //! Path where the populated file (containing the json with all the settings actually used for the simulation)
    //! is going to be saved. Empty string if the file should not be saved.
    path populatedFile = "";
};

//! Create a `json` object from a shared pointer to a `ApplicationOptions` object.
void to_json( json& jsonObject, const boost::shared_ptr< ApplicationOptions >& spiceSettings );

//! Create a shared pointer to a `ApplicationOptions` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< ApplicationOptions >& spiceSettings );

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_OPTIONS_H

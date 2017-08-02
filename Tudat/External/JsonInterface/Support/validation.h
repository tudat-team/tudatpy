/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_VALIDATION_H
#define TUDAT_JSONINTERFACE_VALIDATION_H

#include "valueAccess.h"
#include "valueConversions.h"

namespace tudat
{

namespace json_interface
{

/// ExceptionResponseType

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


/// ValidationSettings

class ValidationSettings
{
public:
    //! Constructor.
    ValidationSettings( ) { }

    //! Response to a "default value used for missing key" event.
    ExceptionResponseType usingDefaultValueForMissingKey = continueSilently;

    //! Response to a "unused key" event.
    ExceptionResponseType unusedKey = printWarning;
};

//! Create a `json` object from a shared pointer to a `ValidationSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< ValidationSettings >& spiceSettings );

//! Create a shared pointer to a `ValidationSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< ValidationSettings >& spiceSettings );

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_VALIDATION_H

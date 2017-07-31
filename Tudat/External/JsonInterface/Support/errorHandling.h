/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_ERRORHANDLING_H
#define TUDAT_JSONINTERFACE_ERRORHANDLING_H

#include <map>

#include <boost/filesystem.hpp>
#include <boost/core/demangle.hpp>

#include "json/src/json.hpp"
using json = nlohmann::json;

#include "keys.h"
#include "utilities.h"

#include <Tudat/InputOutput/basicInputOutput.h>

namespace tudat
{

namespace json_interface
{

//! -DOC
class ValueAccessError : public std::runtime_error
{
public:
    //! Constructor.
    ValueAccessError( const std::string& errorMessage, const KeyPath& keyPath )
        : runtime_error( errorMessage.c_str( ) ), keyPath( keyPath ) { }

    //! Key path.
    const KeyPath keyPath;

    //! Error message.
    virtual const char* what( ) const throw( )
    {
        std::ostringstream stream;
        stream << runtime_error::what( ) << ": " << keyPath;
        std::cerr << stream.str( ).c_str( ) << std::endl;  // FIXME
        return stream.str( ).c_str( );
    }
};

//! -DOC
class UndefinedKeyError : public ValueAccessError
{
public:
    //! Constructor.
    UndefinedKeyError( const KeyPath& keyPath ) : ValueAccessError( "Undefined key", keyPath ) { }
};

//! -DOC
template< typename ExpectedValueType >
class IllegalValueError : public ValueAccessError
{
public:
    //! Constructor.
    IllegalValueError( const KeyPath& keyPath, const json& value )
        : ValueAccessError( "Illegal value for key", keyPath ), value( value ) { }

    //! Associated (illegal) value.
    json value;

    //! Error message.
    virtual const char* what( ) const throw( )
    {
        std::cerr << "Could not convert value to expected type "
                  << boost::core::demangled_name( typeid( ExpectedValueType ) ) << std::endl;
        std::ostringstream stream;
        stream << ValueAccessError::what( ) << " = " << value.dump( 2 );
        // std::cerr << stream.str( ).c_str( ) << std::endl;  // FIXME
        return stream.str( ).c_str( );
    }
};


//! -DOC
class ReportableBugError : public std::exception
{
public:
    //! Constructor.
    ReportableBugError(
            const std::string& errorMessage = "Internal Tudat error. Please, report this bug by using this link: ",
            const std::string& reportLink = "https://github.com/Tudat/tudat/issues" )
        : std::exception( ), errorMessage( errorMessage ), reportLink( reportLink ) { }

    //! Full error message.
    virtual const char* what( ) const throw( )
    {
        std::cerr << errorMessage << std::endl << std::endl << reportLink << std::endl << std::endl;
        std::ostringstream stream;
        stream << std::exception::what( ) << std::endl;
        return stream.str( ).c_str( );
    }

protected:
    //! Error message.
    std::string errorMessage;

    //! Report link.
    std::string reportLink;
};


//! -DOC
class AutoReportableBugError : public ReportableBugError
{
private:
    //! Pre-filled title for new issue.
    std::string issueTitle;

    //! Pre-filled comment for new issue.
    std::string issueBody;

    void updateReportLink( )
    {
        reportLink = "https://github.com/Tudat/tudat/issues/new";
        reportLink += "?title=";
        reportLink += url_encode( issueTitle );
        reportLink += "&body=";
        reportLink += url_encode( issueBody );
    }

public:
    //! Constructor.
    AutoReportableBugError(
            const std::string& errorMessage = "Internal Tudat error. Please, report this bug by using this link: ",
            const std::string& issueTitle = "", const std::string& issueBody = "" )
        : ReportableBugError( errorMessage ), issueTitle( issueTitle ), issueBody( issueBody )
    {
        updateReportLink( );
    }

    void setIssueTitle( const std::string& str )
    {
        issueTitle = str;
        updateReportLink( );
    }

    void setIssueBody( const std::string& str )
    {
        issueBody = str;
        updateReportLink( );
    }

    void setIssue( const std::string& title, const std::string& body )
    {
        issueTitle = title;
        issueBody = body;
        updateReportLink( );
    }
};


//! -DOC
template< typename T >
class NullPointerError : public AutoReportableBugError
{
public:
    //! Constructor
    NullPointerError( ) : AutoReportableBugError( )
    {
        const std::string typeName = boost::core::demangled_name( typeid( T ) );
        errorMessage =
                "Null-pointer of type " + typeName +
                " is not allowed to be NULL.\n"
                "If you think this is not your fault, consider reporting this bug by using this link: ";

        const std::string title = "Null-pointer of type " + typeName;
        const std::string body =
                "I encountered this error when using Tudat's `json_interface`:\n```\n" + errorMessage + "\n```\n";

        setIssue( title, body );
    }
};



/// ENUMS

//! -DOC
class UnknownEnumError : public std::runtime_error
{
public:
    UnknownEnumError( ) : runtime_error( "Unknown conversion between enum and string." ) { }
};

//! -DOC
template< typename EnumType >
EnumType enumFromString( const std::string& stringValue,
                         const std::map< EnumType, std::string >& stringValues )
{
    for ( auto entry : stringValues )
    {
        if ( stringValue == entry.second )
        {
            return entry.first;
        }
    }
    std::cerr << "Unknown string \"" << stringValue << "\" for enum " <<
                 boost::core::demangled_name( typeid( EnumType ) ) << std::endl;
    std::cerr << "Recognized strings:" << std::endl;
    for ( auto entry : stringValues )
    {
        std::cerr << "  " << entry.second << " = " << entry.first << std::endl;
    }
    throw UnknownEnumError( );
}

//! -DOC
template< typename EnumType >
std::string stringFromEnum( const EnumType enumValue, const std::map< EnumType, std::string >& stringValues )
{
    if ( stringValues.count( enumValue ) > 0 )
    {
        return stringValues.at( enumValue );
    }
    else
    {
        std::cerr << "Unknown string representation for enum value "
                  << boost::core::demangled_name( typeid( EnumType ) ) << "::" << enumValue << std::endl;
        throw UnknownEnumError( );
    }
}

//! -DOC
template< typename T >
class UnsupportedEnumError : public AutoReportableBugError
{
public:
    //! Constructor
    UnsupportedEnumError( const T enumValue, const std::map< T, std::string >& stringValues )
        : AutoReportableBugError( )
    {
        const std::string typeName = boost::core::demangled_name( typeid( T ) );
        const std::string stringValue = stringFromEnum( enumValue, stringValues );

        errorMessage =
                "The value \"" + stringValue + "\" for " + typeName +
                "\" is not supported by `json_interface` yet.\n"
                "If you think that support should be added, use this link to request it: ";

        const std::string title = "Add support for \"" + stringValue + "\" in json_interface";
        const std::string body =
                "The value \"" + stringValue + "\" for `" + typeName +
                "` is currently not supported by `json_interface`.\n\n"
                "I think that adding support for it would be a nice addition to Tudat.\n";

        setIssue( title, body );
    }
};


//! -DOC
template< typename T >
class UnimplementedEnumError : public AutoReportableBugError
{
public:
    //! Constructor
    UnimplementedEnumError( const T enumValue, const std::map< T, std::string >& stringValues )
        : AutoReportableBugError( )
    {
        const std::string typeName = boost::core::demangled_name( typeid( T ) );
        const std::string stringValue = stringFromEnum( enumValue, stringValues );

        errorMessage =
                "The value \"" + stringValue + "\" for " + typeName +
                " is marked as supported by `json_interface`, but no implementation was found.\n"
                "This is an internal bug of Tudat. Please, report it by using this link: ";

        const std::string title = "Missing implementation for \"" + stringValue + "\" in json_interface";
        const std::string body =
                "The value \"" + stringValue + "\" for `" + typeName +
                "` is marked as supported by `json_interface`, but no implementation was found.\n";

        setIssue( title, body );
    }
};


//! -DOC
template< typename EnumType >
json handleUnimplementedEnumValueToJson( const EnumType enumValue,
                                         const std::map< EnumType, std::string >& stringValues,
                                         const std::vector< EnumType >& unssupportedValues )
{
    if ( contains( unssupportedValues, enumValue ) )
    {
        UnsupportedEnumError< EnumType > error( enumValue, stringValues );
        error.what( );
    }
    else
    {
        UnimplementedEnumError< EnumType > error( enumValue, stringValues );
        error.what( );
    }
    return json( );
}


//! -DOC
template< typename EnumType >
void handleUnimplementedEnumValueFromJson( const EnumType enumValue,
                                           const std::map< EnumType, std::string >& stringValues,
                                           const std::vector< EnumType >& unssupportedValues )
{
    if ( contains( unssupportedValues, enumValue ) )
    {
        UnsupportedEnumError< EnumType > error( enumValue, stringValues );
        error.what( );
    }
    else
    {
        UnimplementedEnumError< EnumType > error( enumValue, stringValues );
        error.what( );
    }
    throw UnknownEnumError( );
}


template< typename T >
void enforceNonNullPointer( const boost::shared_ptr< T >& pointer )
{
    if ( pointer )
    {
        return;
    }
    else
    {
        throw NullPointerError< T >( );
    }
}

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_ERRORHANDLING_H

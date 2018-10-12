/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include <json/src/json.hpp>

#include "Tudat/JsonInterface/Support/keys.h"
#include "Tudat/JsonInterface/Support/utilities.h"

#include "Tudat/InputOutput/basicInputOutput.h"
namespace tudat
{

namespace json_interface
{

//! Possible responses to an exception during validation phase in Tudat apps that use json_interface.
enum ExceptionResponseType
{
    continueSilently,
    printWarning,
    throwError
};

//! Class for errors generated during the retrieval of a value from a `json` object.
/*!
 * Class for errors generated during the retrieval of a value from a `json` object.
 */
class ValueAccessError : public std::runtime_error
{
public:
    //! Constructor.
    /*!
     * Constructor.
     * \param errorMessage The first part of the message to be printed.
     * \param keyPath Key path trying to access / accessed when the error was generated.
     */
    ValueAccessError( const std::string& errorMessage, const KeyPath& keyPath )
        : runtime_error( errorMessage.c_str( ) ), keyPath( keyPath ) { }

    //! Key path trying to access / accessed when the error was generated.
    KeyPath keyPath;

    //! Full error message.
    virtual const char* what( ) const throw( )
    {
        std::ostringstream stream;
        stream << runtime_error::what( ) << ": " << keyPath;
        std::cerr << stream.str( ).c_str( ) << std::endl;  // FIXME
        return stream.str( ).c_str( );
    }

    //! Whether `this` was generated when trying to access \p keyPath.
    /*!
     * @copybrief wasTriggeredByMissingValueAt
     * \param keyPath The key path to which the instance's key path is to be compared.
     * \return @copybrief wasTriggeredByMissingValueAt
     */
    bool wasTriggeredByMissingValueAt( const KeyPath& keyPath ) const
    {
        return this->keyPath.back( ) == keyPath.back( );
    }

    //! Rethrow `this` if `this` was not generated when trying to access \p keyPath.
    /*!
     * @copybrief rethrowIfNotTriggeredByMissingValueAt
     * \param keyPath The key path to which the instance's key path is to be compared.
     */
    void rethrowIfNotTriggeredByMissingValueAt( const KeyPath& keyPath ) const
    {
        if ( ! wasTriggeredByMissingValueAt( keyPath ) )
        {
            throw *this;
        }
    }
};

//! Class for unrecognized errors generated during the retrieval of a value from a `json` object.
/*!
 * Class for unrecognized errors generated during the retrieval of a value from a `json` object.
 */
class UnrecognizedValueAccessError : public ValueAccessError
{
public:
    //! Constructor.
    /*!
     * Constructor.
     * \param keyPath Key path trying to access when the error was generated.
     * \param expectedTypeInfo Type id info of the expected type.
     */
    UnrecognizedValueAccessError( const KeyPath& keyPath, const std::type_info& expectedTypeInfo ) :
        ValueAccessError( "Unrecognized error when trying to create object of type " +
                          boost::core::demangled_name( expectedTypeInfo ) + " from key", keyPath ) { }
};

//! Class for errors generated when trying to access a key from a `json` object that does not exist.
/*!
 * Class for errors generated when trying to access a key from a `json` object that does not exist.
 */
class UndefinedKeyError : public ValueAccessError
{
public:
    //! Constructor.
    /*!
     * Constructor.
     * \param keyPath Key path trying to access when the error was generated.
     */
    UndefinedKeyError( const KeyPath& keyPath ) : ValueAccessError( "Undefined key", keyPath ) { }

    //! Rethrow `this` if default values are not allowed, or print a message if requested by user.
    /*!
     * @copybrief handleUseOfDefaultValue
     * \param defaultValue `json` representation of the dafault value to be used.
     * \param response The response to usage of dafault values.
     */
    void handleUseOfDefaultValue( const nlohmann::json& defaultValue, const ExceptionResponseType& response ) const
    {
        if ( ! containsAnyOf( keyPath, SpecialKeys::objectContaining ) )
        {
            if ( response != continueSilently )
            {
                if ( response == throwError )
                {
                    throw *this;
                }
                else
                {
                    std::cerr << "Using default value for key: " << keyPath << " = " << defaultValue << std::endl;
                }
            }
        }
    }

    //! Rethrow `this` if NaN default values are not allowed and the provided default value is NaN.
    /*!
     * @copybrief rethrowIfNaNNotAllowed
     * \param allowNaN Whether NaN default values are allowed.
     * \param defaultValue The default value to be checked.
     */
    template< typename NumberType >
    void rethrowIfNaNNotAllowed( bool allowNaN, const NumberType& defaultValue ) const
    {
        if ( ! allowNaN && isNaN( defaultValue ) )
        {
            throw *this;
        }
    }
};

//! Class for errors generated when trying to convert a `json` object to another type.
/*!
 * Class for errors generated when trying to convert a `json` object to another type.
 */
class IllegalValueError : public ValueAccessError
{
public:
    //! Constructor.
    /*!
     * Constructor.
     * \param keyPath Key path \p value was retrieved from.
     * \param value The `json` object that was going to be converted to the expected type.
     * \param expectedTypeInfo Type id info of the expected type.
     */
    IllegalValueError( const KeyPath& keyPath, const nlohmann::json& value, const std::type_info& expectedTypeInfo ) :
        ValueAccessError( "Illegal value for key", keyPath ), value( value ),
        expectedTypeName( boost::core::demangled_name( expectedTypeInfo ) ) { }

    //! Associated (illegal) value.
    nlohmann::json value;

    //! Expected type name.
    std::string expectedTypeName;

    //! Full error message.
    virtual const char* what( ) const throw( )
    {
        std::cerr << "Could not convert value to expected type " << expectedTypeName << std::endl;
        std::ostringstream stream;
        stream << ValueAccessError::what( );
        return stream.str( ).c_str( );
    }
};

//! Class for errors that print a report bug link.
/*!
 * Class for errors that print a report bug link.
 */
class ReportableBugError : public std::exception
{
public:
    //! Constructor.
    /*!
     * Constructor.
     * \param errorMessage @copydoc ReportableBugError::errorMessage
     * \param reportURL @copydoc ReportableBugError::reportURL
     */
    ReportableBugError(
            const std::string& errorMessage = "Internal Tudat error. Please, report this bug by using this link: ",
            const std::string& reportURL = "https://github.com/Tudat/tudat/issues" )
        : std::exception( ), errorMessage( errorMessage ), reportURL( reportURL ) { }

    //! Full error message.
    virtual const char* what( ) const throw( )
    {
        std::cerr << errorMessage << std::endl << std::endl << reportURL << std::endl << std::endl;
        std::ostringstream stream;
        stream << std::exception::what( ) << std::endl;
        return stream.str( ).c_str( );
    }

protected:
    //! The message to print before the error link.
    std::string errorMessage;

    //! The URL that can be used to report the error.
    std::string reportURL;
};


//! Class for errors that print a report bug link to open an issue on GitHub (optionally pre-filled).
/*!
 * Class for errors that print a report bug link to open an issue on GitHub (optionally pre-filled).
 */
class AutoReportableBugError : public ReportableBugError
{
public:
    //! Constructor.
    /*!
     * Constructor.
     * \param errorMessage First part of the message.
     * \param issueTitle @copydoc AutoReportableBugError::issueTitle
     * \param issueBody @copydoc AutoReportableBugError::issueBody
     */
    AutoReportableBugError(
            const std::string& errorMessage = "Internal Tudat error. Please, report this bug by using this link: ",
            const std::string& issueTitle = "", const std::string& issueBody = "" )
        : ReportableBugError( errorMessage ), issueTitle( issueTitle ), issueBody( issueBody )
    {
        updateReportURL( );
    }

    //! Change the issue's title and update the report URL.
    /*!
     * @copybrief setIssueTitle
     * \param title New issue's title.
     */
    void setIssueTitle( const std::string& title )
    {
        issueTitle = title;
        updateReportURL( );
    }

    //! Change the issue's body and update the report URL.
    /*!
     * @copybrief setIssueBody
     * \param body New issue's body.
     */
    void setIssueBody( const std::string& body )
    {
        issueBody = body;
        updateReportURL( );
    }

    //! Change the issue's title and body and update the report URL.
    /*!
     * @copybrief setIssue
     * \param title New issue's title.
     * \param body New issue's body.
     */
    void setIssue( const std::string& title, const std::string& body )
    {
        issueTitle = title;
        issueBody = body;
        updateReportURL( );
    }

private:
    //! Pre-filled title for the new issue to be opened on GitHub.
    std::string issueTitle;

    //! Pre-filled title for the new issue to be opened on GitHub.
    std::string issueBody;

    //! Update the report URL.
    /*!
     * Update the report URL for the current \p issueTitle and \p issueBody.
     */
    void updateReportURL( )
    {
        reportURL = "https://github.com/Tudat/tudat/issues/new";
        reportURL += "?title=";
        reportURL += url_encode( issueTitle );
        reportURL += "&body=";
        reportURL += url_encode( issueBody );
    }
};


// nullptr POINTERS

//! Class for errors generated when a pointer that should not be `nullptr` is `nullptr`.
/*!
 * Class for errors generated when a pointer to an object of type `T` that should not be `nullptr` is `nullptr`.
 */
template< typename T >
class nullptrPointerError : public AutoReportableBugError
{
public:
    //! Constructor
    /*!
     * Empty constructor.
     */
    nullptrPointerError( ) : AutoReportableBugError( )
    {
        const std::string typeName = boost::core::demangled_name( typeid( T ) );
        errorMessage =
                "nullptr-pointer of type " + typeName +
                " is not allowed to be nullptr.\n"
                "If you think this is not your fault, consider reporting this bug by using this link: ";

        const std::string title = "nullptr-pointer of type " + typeName;
        const std::string body =
                "I encountered this error when using Tudat's `json_interface`:\n```\n" + errorMessage + "\n```\n";

        setIssue( title, body );
    }
};

//! Check that a pointer is not `nullptr`.
/*!
 * Check that a pointer to an object of type `T` is not `nullptr`.
 * \param pointer The pointer that is not allowed to be `nullptr`.
 * \throws nullptrPointerError< T > If \p pointer is `nullptr`.
 */
template< typename T >
void assertNonnullptrPointer( const std::shared_ptr< T >& pointer )
{
    if ( pointer )
    {
        return;
    }
    else
    {
        throw nullptrPointerError< T >( );
    }
}



// ENUMS

//! Class for errors generated when trying to convert between `enum` and `std::string`.
/*!
 * Class for errors generated when trying to convert between `enum` and `std::string`.
 */
class UnknownEnumError : public std::runtime_error
{
public:
    //! Constructor.
    /*!
     * Empty constructor.
     */
    UnknownEnumError( ) : runtime_error( "Unknown conversion between enum and string." ) { }
};


// enum <-> std::string

//! Get an `enum` value of type `EnumType` from a `std::string`.
/*!
 * @copybrief \enumFromString
 * \param stringValue The string to be converted to `EnumType`.
 * \param stringValues Map containing the string representation for the values of `EnumType`.
 * \return The `EnumType` corresponding to \p stringValue.
 * \throws UnknownEnumError If no entry with value \p stringValue is found in \p stringValues.
 */
template< typename EnumType >
EnumType enumFromString( const std::string& stringValue, const std::map< EnumType, std::string >& stringValues )
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
        std::cerr << "  " << entry.second << std::endl;
    }
    throw UnknownEnumError( );
}

//! Get a `std::string` from an `enum` value of type `EnumType`.
/*!
 * @copybrief \stringFromEnum
 * \param enumValue The `EnumType` to be converted to string.
 * \param stringValues Map containing the string representation for the values of `EnumType`.
 * \return The `std::string` corresponding to \p enumValue.
 * \throws UnknownEnumError If no entry with key \p enumValue is found in \p stringValues.
 */
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


//! Class for errors generated when trying to use a value for an `enum` of type `T` that is not supported by
//! json_interface.
/*!
 * Class for errors generated when trying to use a value for an `enum` of type `T` that is not supported by
 * json_interface.
 */
template< typename T >
class UnsupportedEnumError : public std::runtime_error
{
public:
    //! Constructor.
    /*!
     * Constructor.
     * \param enumValue The `enum` value of type `T` that is not supported.
     * \param stringValues Map containing the string representation for the values of `EnumType`
     * (including \p enumValue).
     */
    UnsupportedEnumError( const T enumValue, const std::map< T, std::string >& stringValues )
        : std::runtime_error( "The value \"" + stringFromEnum( enumValue, stringValues ) + "\" for " +
                              boost::core::demangled_name( typeid( T ) ) +
                              "\" is not supported directly by the `json_interface`.\n" +
                              "Write your own JSON-based C++ application if you want to use this Tudat feature " +
                              "in combination with JSON input files." ) { }
};

//! Class for errors generated when trying to use a value for an `enum` of type `EnumType` that is marked as supported
//! by json_interface but for which no implementation was found.
/*!
 * Class for errors generated when trying to use a value for an `enum` of type `EnumType` that is marked as supported
 * by json_interface but for which no implementation was found.
 */
template< typename EnumType >
class UnimplementedEnumError : public AutoReportableBugError
{
public:
    //! Constructor.
    /*!
     * Constructor.
     * \param enumValue The `enum` value of type `EnumType` that is not supported.
     * \param stringValues Map containing the string representation for the values of `EnumType`
     * (including \p enumValue).
     */
    UnimplementedEnumError( const EnumType enumValue, const std::map< EnumType, std::string >& stringValues )
        : AutoReportableBugError( )
    {
        const std::string typeName = boost::core::demangled_name( typeid( EnumType ) );
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

//! Handle an unimplemented `enum` value of type `EnumType`.
/*!
 * @copybrief handleUnimplementedEnumValue
 * \param enumValue The `enum` value of type `T` that is not implemented.
 * \param stringValues Map containing the string representation for the values of `EnumType`
 * (including \p enumValue).
 * \param unssupportedValues Vector of values of `EnumType` that are not supported by json_interface.
 * \throws UnsupportedEnumError<EnumType> If \p enumValue is contained in \p unssupportedValues.
 * \throws UnimplementedEnumError<EnumType> If \p enumValue is not contained in \p unssupportedValues.
 */
template< typename EnumType >
void handleUnimplementedEnumValue( const EnumType enumValue,
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
}

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_ERRORHANDLING_H

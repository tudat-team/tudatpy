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

#include <boost/make_shared.hpp>

#include "Tudat/InputOutput/fieldValue.h"

namespace tudat
{
namespace input_output
{

//! Create a FieldValue containing type, string content and transformation of field.
FieldValue::FieldValue( const FieldType& fieldType, const std::string& fieldContent,
                        const boost::shared_ptr< FieldTransform > transformer )
    : type( fieldType ), rawField( fieldContent ), transform ( transformer )
{ }

//! Get transformed field content.
const std::string& FieldValue::getTransformed( )
{
    loadTransformation( );
    return transformedField;
}

//! Get raw field content.
const std::string& FieldValue::getRaw( ) { return rawField; }

} // namespace input_output
} // namespace tudat

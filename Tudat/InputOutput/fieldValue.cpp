/*   Copyright (c) 2010-2011 Delft University of Technology.
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
 *      120326    D. Dirkx          Code checked, minor layout changes.
 */

#include "Tudat/InputOutput/fieldValue.h"

namespace tudat
{
namespace input_output
{

//! Create a FieldValue containing type, string content and transformation of field.
FieldValue::FieldValue( FieldType& fieldType, std::string& fieldContent,
                        boost::shared_ptr<FieldTransform> transformer ) :
    // Set the field type
    type( fieldType ),
    // Copy the raw field into the object
    rawField( fieldContent ),
    // Set the corresponding transformer
    transform (transformer ) { }

//! Get value of field content in SI units.
boost::shared_ptr<std::string> FieldValue::get( )
{
    // Check if field has a unit transformation.
    if ( transform.get( ) )
        return  transform->transform( getRaw( ) );
    else
        return  getRaw( );
}

//! Get raw field content.
boost::shared_ptr<std::string> FieldValue::getRaw( )
{
    return  boost::shared_ptr<std::string>( new std::string( rawField ) );
}

//! Operator to get value of field content in SI units.
boost::shared_ptr<std::string> FieldValue::operator( ) ( )
{
    return get( );
}

} // namespace input_output
} // namespace tudat

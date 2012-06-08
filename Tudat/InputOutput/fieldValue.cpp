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
 *      120326    D. Dirkx          Code checked, minor layout changes.
 *
 *    References
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
    : type( fieldType ), // Set the field type.
      rawField( fieldContent ), // Copy the raw field into the object.
      transform ( transformer )  // Set the corresponding transformer.
{ }

//! Get value of field content in SI units.
boost::shared_ptr< std::string > FieldValue::get( )
{
    // Check if field has a unit transformation.
    if ( transform.get( ) )
    {
        return  transform->transform( getRaw( ) );
    }
    else
    {
        return getRaw( );
    }
}

//! Get raw field content.
boost::shared_ptr< std::string > FieldValue::getRaw( )
{
    return boost::make_shared< std::string >( rawField );
}

//! Operator to get value of field content in SI units.
boost::shared_ptr< std::string > FieldValue::operator( ) ( )
{
    return get( );
}

} // namespace input_output
} // namespace tudat

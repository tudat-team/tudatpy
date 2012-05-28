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
 *      111103    S. Billemont      First creation of code.
 *      120217    D.J. Gondelach    Code check.
 *      120326    D. Dirkx          Code checked, minor layout changes.
 *
 *    References
 */

#ifndef TUDAT_FIELD_VALUE_H
#define TUDAT_FIELD_VALUE_H

#include <string>

#include <boost/shared_ptr.hpp>

#include "Tudat/InputOutput/fieldTransform.h"
#include "Tudat/InputOutput/fieldType.h"

namespace tudat
{
namespace input_output
{

//! Field value class.
class FieldValue
{
public:
    //! Create a FieldValue containing type, string content and transformation of field.
    /*! Create a FieldValue which contains the type, the content as a string
     *  and unit transformation of the field.
     */
    FieldValue( FieldType& type, std::string& field,
                boost::shared_ptr< FieldTransform > transformer
                = boost::shared_ptr< FieldTransform >( ) );

    //! FieldType of the FieldValue.
    FieldType type;

    //! Get value of field content in SI units.
    boost::shared_ptr< std::string > get( );

    //! Get raw field content.
    boost::shared_ptr< std::string > getRaw( );

    //! Operator to get value of field content in SI units.
    boost::shared_ptr< std::string > operator( ) ( );

protected:

private:

    //! String of field value in raw data format.
    std::string rawField;

    //! Pointer to unit transformation equation.
    boost::shared_ptr< FieldTransform > transform;
};

} // namespace input_output
} // namespace tudat

#endif // TUDAT_FIELD_VALUE_H

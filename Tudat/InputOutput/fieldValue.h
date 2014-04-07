/*    Copyright (c) 2010-2014, Delft University of Technology
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
 *      120217    D.J. Gondelach    Code check.
 *      120326    D. Dirkx          Code checked, minor layout changes.
 *      130121    K. Kumar          Added shared-ptr typedef.
 *      130301    S. Billemont      Update FieldValue definition.
 *      131221    K. Kumar          Fixed Doxygen comments.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_FIELD_VALUE_H
#define TUDAT_FIELD_VALUE_H

#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/smart_ptr/make_shared.hpp>

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
    /*!
     * Create a FieldValue which contains the type, the content as a string and unit transformation
     * of the field.
     * \param type Field type.
     * \param field Field data.
     * \param transformer Transformation to apply to field data.
     */
    FieldValue( const FieldType& type, const std::string& field,
                const boost::shared_ptr< FieldTransform > transformer
                = boost::shared_ptr< FieldTransform >( ) );

    //! FieldType of the FieldValue.
    const FieldType type;

    //! Get value of field content in SI units.
    template< typename T >
    T get( )
    {
        loadTransformation( );
        return boost::lexical_cast< T >( transformedField );
    }

    //! Get a shared pointer to the value of field content in SI units.
    template< typename T >
    boost::shared_ptr< T > getPointer( )
    {
        loadTransformation( );
        return boost::make_shared< T >( boost::lexical_cast<T>( transformedField ) );
    }

    //! Get transformed field content.
    const std::string& getTransformed( );

    //! Get raw field content.
    const std::string& getRaw( );

protected:

    //! Load the transformed string if required.
    inline void loadTransformation( )
    {
        // Check if transformedField exists already (lazy cache value).
        if ( transformedField.empty( ) )
        {
            // Doesn't exist yet, load the transform.
            transformedField = transform.get( ) ? 
                *transform->transform( rawField ).get( ) : rawField;
        }
    }

private:

    //! String of field value in raw data format.
    const std::string   rawField;

    //! String of field value in the transformed data format.
    std::string         transformedField;

    //! Pointer to unit transformation equation.
    FieldTransformPointer transform;
};

//! Typedef for shared-pointer to FieldValue object.
typedef boost::shared_ptr< FieldValue > FieldValuePointer;

} // namespace input_output
} // namespace tudat

#endif // TUDAT_FIELD_VALUE_H

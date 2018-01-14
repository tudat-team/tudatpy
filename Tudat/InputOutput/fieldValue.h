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

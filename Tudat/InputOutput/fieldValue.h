/*   Copyright (c) 2010-2012 Delft University of Technology.
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

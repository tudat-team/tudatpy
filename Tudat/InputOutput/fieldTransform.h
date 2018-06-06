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

#ifndef TUDAT_FIELD_TRANSFORM_H
#define TUDAT_FIELD_TRANSFORM_H

#include <memory>

#include <string>

namespace tudat
{
namespace input_output
{

//! Field transform base class.
/*!
 * This abstract class belongs to the parser-extractor architecture implemented in Tudat. This base
 * class can be used to derive specific transformation classes that take strings and return
 * shared-pointers to transformed strings.
 * \sa Extractor, Parser
 */
class FieldTransform
{
public:

    //! Default destructor.
    virtual ~FieldTransform( ) { }

    //! Transform input string.
    /*! 
     * Returns a transformed string.
     * \param input Input string.
     * \return Shared-pointer to transformed string.
     */
    virtual std::shared_ptr< std::string > transform( const std::string& input ) = 0;

protected:

private:
};

//! Typedef for shared-pointer to FieldTransform object.
typedef std::shared_ptr< FieldTransform > FieldTransformPointer;

} // namespace input_output
} // namespace tudat

#endif // TUDAT_FIELD_TRANSFORM_H

/*    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      120326    D. Dirkx          Code checked, minor layout changes.
 */

#ifndef FIELDTRANSFORM_H
#define FIELDTRANSFORM_H

#include <boost/shared_ptr.hpp>
#include <string>

namespace tudat
{
namespace input_output
{

class FieldTransform
{

public:

    //! Default constructor.
    FieldTransform( ) { }

    //! Default destructor.
    virtual ~FieldTransform( ) { }

    //! Transform input string.
    /*!
     * Returns a transformed string.
     *
     * \param input Input string.
     */
    virtual boost::shared_ptr<std::string> transform( boost::shared_ptr<std::string> input ) = 0;

};
} // namespace input_output
} // namespace tudat
#endif

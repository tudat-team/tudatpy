/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      120326    D. Dirkx          Code checked, minor layout changes.
 *      130121    K. Kumar          Added shared-ptr typedef.
 *      131221    K. Kumar          Fixed missing Doxygen comments.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_FIELD_TRANSFORM_H
#define TUDAT_FIELD_TRANSFORM_H

#include <boost/shared_ptr.hpp>

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
    virtual boost::shared_ptr< std::string > transform( const std::string& input ) = 0;

protected:

private:
};

//! Typedef for shared-pointer to FieldTransform object.
typedef boost::shared_ptr< FieldTransform > FieldTransformPointer;

} // namespace input_output
} // namespace tudat

#endif // TUDAT_FIELD_TRANSFORM_H

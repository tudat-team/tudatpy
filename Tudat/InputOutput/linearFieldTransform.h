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
 *      120906    T. Secretin       Reviewed code. Moved transform implementation to .cpp file.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_LINEAR_FIELD_TRANSFORM_H
#define TUDAT_LINEAR_FIELD_TRANSFORM_H

#include <string>

#include <boost/shared_ptr.hpp>

#include "Tudat/InputOutput/fieldTransform.h"

namespace tudat
{
namespace input_output
{

class LinearFieldTransform : public FieldTransform
{
public:

    //! Constructor of the linear field transform.
    /*!
     * Constructor of the linear field transformation. Linear transformation is of the form:
     * \f$y=a*x+b\f$, where a is the slope and b is the intercept.
     * \param a Slope of the linear field transform.
     * \param b Intercept of the linear field transform.
     */
    LinearFieldTransform( const double aSlope, const double anIntercept )
        : slope( aSlope ), intercept( anIntercept )
    { }

    //! Default destructor.
    ~LinearFieldTransform( ) { }

    //! Transform input string.
    /*!
     * Returns a transformed string, according to the linear transformation:
     * result = slope * input + intercept.
     * \param input Input string.
     * \return Transformed string.
     */
    boost::shared_ptr< std::string > transform( boost::shared_ptr< std::string > input );

protected:

private:

    //! Slope of the transformation.
    const double slope;

    //! Intercept of the transformation.
    const double intercept;
};

} // namespace input_output
} // namespace tudat

#endif // TUDAT_LINEAR_FIELD_TRANSFORM_H

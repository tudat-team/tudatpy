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
 *      120716    D. Dirkx          Creation of file.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_INTERPOLATOR_H
#define TUDAT_INTERPOLATOR_H

#include <vector>

namespace tudat
{
namespace interpolators
{

//! Base class for interpolator.
/*!
 * Base class for the interpolators included in Tudat, the dependent and independent variable
 * types are specified as user parameters, as are the number of dimensions. The number of
 * dimensions is chosen as a template parameter, so that a boost multi_array of this
 * dimension can be used as a member variable of derived classes.
 * \tparam IndependentVariableType Type of independent variable(s).
 * \tparam IndependentVariableType Type of dependent variable.
 * \tparam numberOfDimensions Number of independent directions for independent variables.
 */
template< typename IndependentVariableType, typename DependentVariableType,
          int numberOfDimensions >
class Interpolator
{
public:

    //! Destructor.
    /*!
     * Destructor.
     */
    virtual ~Interpolator( ) { }

    //! Interpolate.
    /*!
     * This function performs the interpolation. It must be implemented in derived classes.
     * \param independentVariableValues Vector of values of independent variables at which
     *          the value of the dependent variable is to be determined.
     * \param Interpolated value of dependent variable.
     */
    virtual DependentVariableType interpolate( const std::vector< IndependentVariableType >&
                                               independentVariableValues ) = 0;

};

} // namespace interpolators
} // namespace tudat

#endif // TUDAT_INTERPOLATOR_H

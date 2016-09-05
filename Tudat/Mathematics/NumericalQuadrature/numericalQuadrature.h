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
 *      120716    D. Dirkx          Creation of file.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_NUMERICAL_QUADRATURE_H
#define TUDAT_NUMERICAL_QUADRATURE_H

#include <vector>
#include <Eigen/Core>

namespace tudat
{
namespace numerical_quadrature
{

//! Base class for numerical quadrature.
/*!
 * Base class for the numerical quadrature methods included in Tudat, the dependent and independent variable
 * types are specified as user parameters, as are the number of dimensions. The number of
 * dimensions is chosen as a template parameter, so that a boost multi_array of this
 * dimension can be used as a member variable of derived classes.
 * \tparam IndependentVariableType Type of independent variable(s).
 * \tparam IndependentVariableType Type of dependent variable.
 */
template< typename IndependentVariableType, typename DependentVariableType >
class NumericalQuadrature
{
public:

    //! Destructor.
    /*!
     * Destructor.
     */
    virtual ~NumericalQuadrature( ) { }

    //! Integrate.
    /*!
     * This function performs the integration. It must be implemented in derived classes.
     * \return Integral of the data.
     */
    virtual DependentVariableType integrate(  ) = 0;

protected:

private:


};

} // namespace numerical_quadrature
} // namespace tudat

#endif // TUDAT_NUMERICAL_QUADRATURE_H

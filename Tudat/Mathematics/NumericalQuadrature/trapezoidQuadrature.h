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

#ifndef TUDAT_TRAPEZOIDAL_INTEGRATOR_H
#define TUDAT_TRAPEZOIDAL_INTEGRATOR_H

#include <vector>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <Tudat/Mathematics/NumericalQuadrature/numericalQuadrature.h>

namespace tudat
{
namespace numerical_quadrature
{


//! Integrate.
/*!
 * This function performs the integration.
 * Source: Burden & Faires
 * \return Integral of the data.
 */
template< typename IndependentVariableType, typename DependentVariableType >
DependentVariableType performTrapezoidalQuadrature(
        const std::vector< IndependentVariableType >& independentVariables,
                                           const std::vector< DependentVariableType >& dependentVariables )
{
    DependentVariableType integral = dependentVariables.at( 0 ) - dependentVariables.at( 0 );
    DependentVariableType h;
    for( unsigned int i = 0 ; i < independentVariables.size() - 1 ; i++ )
    {
        h = independentVariables[i+1] - independentVariables[i];
        integral += h * ( dependentVariables[i+1] + dependentVariables[ i ] ) / ( 2.0 ) ;
    }
    return integral;
}

//! Trapezoid numerical integrator.
/*!
 * Numerical integrator that uses the trapezoid method to compute definite integrals of a dataset.
 * \tparam VariableType Type of independent and dependent variable(s).
 */
template< typename IndependentVariableType, typename DependentVariableType >
class TrapezoidNumericalQuadrature : public NumericalQuadrature< IndependentVariableType , DependentVariableType >
{
public:

    //! Constructor
    TrapezoidNumericalQuadrature( const std::vector< IndependentVariableType >& independentVariables,
                                  const std::vector< DependentVariableType >& dependentVariables)
    {
        independentVariables_ = independentVariables;
        dependentVariables_ = dependentVariables;
    }

    //! Constructor
    TrapezoidNumericalQuadrature( ){ }

    //! Integrate.
    /*!
     * This function performs the integration.
     * Source: Burden & Faires
     * \return Integral of the data.
     */
    DependentVariableType integrate( )
    {
        return performTrapezoidalQuadrature( independentVariables_ , dependentVariables_ );
    }

protected:

private:

};

typedef boost::shared_ptr< TrapezoidNumericalQuadrature< double, double > > TrapezoidNumericalIntegratorPointerd;

} // namespace numerical_quadrature
} // namespace tudat

#endif // TUDAT_TRAPEZOIDAL_INTEGRATOR_H

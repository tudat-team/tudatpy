/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      Numerical differentiation has been implemented using only central difference. In future,
 *      this will be updated to the functionality present in numericalDerivative.h. Numerical
 *      integration for the computeDefiniteIntegral has not been implemented yet. This will be
 *      implemented in future with the Runge-Kutta 4 are the default numerical integrator.
 *
 */

#ifndef TUDAT_BASIC_FUNCTION_H
#define TUDAT_BASIC_FUNCTION_H

#include <limits>
#include <stdexcept>

#include <boost/bind.hpp>
#include <memory>

//#include "Tudat/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h"

#include "Tudat/Mathematics/BasicMathematics/function.h"

namespace tudat
{
namespace basic_mathematics
{

//! Function with basic implementations for derivatives and integrals.
/*!
 * This class fills in the computeDerivative and computeDefiniteIntegral functionality using
 * numerical differentiation and integration. These are slower than their analytical counterparts,
 * but useful when a closed form expression is available.
 */
template< typename IndependentVariable = double, typename DependentVariable = double >
class BasicFunction : public Function< IndependentVariable, DependentVariable >
{
public:

    //! Number of steps to take when performing numerical integration.
    unsigned int integrationSteps;

    //! Constructor with definition of the number of steps.
    BasicFunction( unsigned int numberOfIntegrationSteps = 1000 )
        : integrationSteps( numberOfIntegrationSteps )
    { }

    //! Default destructor.
    virtual ~BasicFunction( ) { }

    //! Derivative (of a given order) of a function.
    /*!
     * Generic calculation of the derivative of a function using numerical techniques.
     * \see Function::derivative
     * TODO: Use the external derivative code, this function is a basic/temporary implementation.
     * for now this is the implementation:
     * \param order The Order of the derivative
     * \param independentVariable The independent variable with respect to which the derivative is taken.
     * \return The numerical derivative of the function.
     */
    virtual DependentVariable computeDerivative( const unsigned int order,
                                                 const IndependentVariable independentVariable )
    {
        if ( order == 0 )
        {
            // Zero th order derivative is the function itself.
            return this->evaluate( independentVariable );
        }

        // Define step size to the left and right.
        double stepSize = sqrt_epsilon_double * independentVariable;

        // Return derivative according to formula:
        // f'(x) = [f(x+h)-f(x-h)]/2h
        return ( computeDerivative( order - 1, independentVariable + stepSize )
                 - computeDerivative( order - 1, independentVariable - stepSize ) )
                / ( 2.0 * stepSize );
    }

    //! Alias for computeDefiniteIntegral, but with an extra variable. Should not be used!
    inline DependentVariable computeDefiniteIntegralMock(
            const unsigned int order, const IndependentVariable lowerBound,
            const IndependentVariable upperbound, const DependentVariable placeHolder )
    {
        return computeDefiniteIntegral( order, lowerBound, upperbound );
    }

    //! Integral (of a given order) of a function from a given lower bound to upper bound
    //! (without constants).
    /*!
     * Generic calculation of the function integral using numerical techniques.
     * TODO: Fix this code, as it doesnt work.
     * \param order Order of integral (i.e. integrate n times, with n = order )
     * \param lowerBound Lower integration bound
     * \param upperbound Upper integration bound
     * \return Nothing as yet; function not yet implemented
     * \see Function::computeDefiniteIntegral.
     */
    virtual DependentVariable computeDefiniteIntegral( const unsigned int order,
                                                       const IndependentVariable lowerBound,
                                                       const IndependentVariable upperbound )
    {
        throw std::runtime_error( "Numerical integration not yet implemented!" );
    }

protected:

    //! Square root of double precision, used to scale the step size.
    /*!
     * The reason not to scale with epsilon itself is to limit the error due to floating point 
     * arithmetic with super small numbers: = sqrt(std::numeric_limits<double>::epsilon()).
     * TODO: Remove when using external derivative code.
     */
    static const double sqrt_epsilon_double;

private:

    //! typedef for a Function that relates a Dependent Variable to an IndependentVariable.
    typedef Function< IndependentVariable, DependentVariable >          Parent;
    //! typedef for a Runge-Kutta 4 integration scheme.
//    typedef numerical_integrators::RungeKutta4Integrator<
//    IndependentVariable, DependentVariable, DependentVariable >         Integrator;
};

// TODO: Remove when using external derivative code.
template< typename IndependentVariable, typename DependentVariable >
const double BasicFunction< IndependentVariable, DependentVariable >::sqrt_epsilon_double
= std::numeric_limits< double >::epsilon( );

// Some handy typedefs.
typedef std::shared_ptr< BasicFunction< double, double > > BasicFunctionPointer;

} // namespace basic_mathematics
} // namespace tudat

#endif // TUDAT_BASIC_FUNCTION_H

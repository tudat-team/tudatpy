/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *      Press, W., Flannery, B., Teukolsky, S., and Vetterling, W., Numerical Recipes in Fortran 77:
 *          The Art of Scientific Computing, 2nd ed. Cambridge University Press, 1992, vol. 1.
 */

#ifndef TUDAT_TRAPEZOIDAL_INTEGRATOR_H
#define TUDAT_TRAPEZOIDAL_INTEGRATOR_H

#include <vector>

#include <memory>

#include <Eigen/Core>

#include "Tudat/Mathematics/NumericalQuadrature/numericalQuadrature.h"
#include "Tudat/Basics/identityElements.h"

namespace tudat
{

namespace numerical_quadrature
{

//! Function to perform numerical quadrature using the trapezoidal method.
/*!
 *  Function to perform numerical quadrature using the trapezoidal method.
 *  \param independentVariables Values of independent variables at which dependentVariables are given.
 *  \param dependentVariables Values of function for which the numerical quadrature is to be computed, given at
 *      independentVariables.
 *  \return Numerical quadrature (integral) of the data provided as input.
 */
template< typename IndependentVariableType, typename DependentVariableType >
DependentVariableType performTrapezoidalQuadrature(
        const std::vector< IndependentVariableType >& independentVariables,
        const std::vector< DependentVariableType >& dependentVariables )
{
    DependentVariableType integral = dependentVariables.at( 0 ) - dependentVariables.at( 0 );
    IndependentVariableType timeStep;
    for( unsigned int i = 0 ; i < independentVariables.size( ) - 1 ; i++ )
    {
        timeStep = independentVariables[ i + 1 ] - independentVariables[ i ];
        integral += timeStep * ( dependentVariables[ i + 1 ] + dependentVariables[ i ] ) / 2.0 ;
    }
    return integral;
}

//! Function to perform numerical quadrature using the extended Simpson's method.
/*!
 *  Function to perform numerical quadrature using the extended Simpson's method (Press et al., 1992). Note
 *  that the spacing of the (in)dependent needs to be constant.
 *  \param constantIndependentVariableStep Constant independent variable step size.
 *  \param dependentVariables Values of function for which the numerical quadrature is to be computed, given at
 *      independentVariables.
 *  \return Numerical quadrature (integral) of the data provided as input.
 */
template< typename IndependentVariableType, typename DependentVariableType >
DependentVariableType performExtendedSimpsonsQuadrature(
        const IndependentVariableType constantIndependentVariableStep,
        const std::vector< DependentVariableType >& dependentVariables )
{
    // Get initial variables
    DependentVariableType integral = dependentVariables.at( 0 ) - dependentVariables.at( 0 );
    unsigned int numberOfVariables = dependentVariables.size( );

    // Check that there are enough elements in the vector
    if ( numberOfVariables > 6 )
    {
        // Create vector of weights
        std::vector< IndependentVariableType > vectorOfWeights = std::vector< IndependentVariableType >( numberOfVariables,
                                                                                                         constantIndependentVariableStep );

        // Add weights
        vectorOfWeights.at( 0 ) *= 3.0 / 8.0;
        vectorOfWeights.at( numberOfVariables - 1 ) *= 3.0 / 8.0;

        vectorOfWeights.at( 1 ) *= 7.0 / 6.0;
        vectorOfWeights.at( numberOfVariables - 2 ) *= 7.0 / 6.0;

        vectorOfWeights.at( 2 ) *= 23.0 / 24.0;
        vectorOfWeights.at( numberOfVariables - 3 ) *= 23.0 / 24.0;

        // Loop over each time step and return result
        for( unsigned int i = 0 ; i < numberOfVariables; i++ )
        {
            integral += vectorOfWeights[ i ] * dependentVariables[ i ];
        }
    }
    else if ( numberOfVariables != 1 )
    {
        // Perform trapezoidal integration instead
        for( unsigned int i = 0 ; i < ( numberOfVariables - 1 ); i++ )
        {
            integral += constantIndependentVariableStep * (
                        dependentVariables[ i + 1 ] + dependentVariables[ i ] ) / 2.0 ;
        }
    }
    else
    {
        integral = IdentityElement::getAdditionIdentity< DependentVariableType >( );
    }
    return integral;
}

//! Trapezoid numerical quadrature wrapper class.
/*!
 *  Numerical method that uses the trapezoid method to compute definite integrals of a dataset.
 */
template< typename IndependentVariableType, typename DependentVariableType >
class TrapezoidNumericalQuadrature : public NumericalQuadrature< IndependentVariableType , DependentVariableType >
{
public:

    //! Constructor.
    /*!
     *  Constructor
     *  \param independentVariables Values of independent variables at which dependentVariables are given
     *  \param dependentVariables Values of function for which the numerical quadrature is to be computed, given at
     *      independentVariables.
     */
    TrapezoidNumericalQuadrature( const std::vector< IndependentVariableType >& independentVariables,
                                  const std::vector< DependentVariableType >& dependentVariables )
    {
        independentVariables_ = independentVariables;
        dependentVariables_ = dependentVariables;
        performQuadrature( );
    }

    //! Function to reset the (in)dependent variable values.
    /*!
     *  Function to reset the (in)dependent variable values.
     *  \param independentVariables Values of independent variables at which dependentVariables are given
     *  \param dependentVariables Values of function for which the numerical quadrature is to be computed, given at
     *      independentVariables.
     */
    void resetData( const std::vector< IndependentVariableType >& independentVariables,
                    const std::vector< DependentVariableType >& dependentVariables)
    {
        independentVariables_ = independentVariables;
        dependentVariables_ = dependentVariables;
        performQuadrature( );
    }

    //! Function to return computed value of the quadrature.
    /*!
     *  Function to return computed value of the quadrature, as computed by last call to performQuadrature.
     *  \return Function to return computed value of the quadrature, as computed by last call to performQuadrature.
     */
    DependentVariableType getQuadrature( )
    {
        return quadratureResult_;
    }

protected:

    //! Function that is called to perform the numerical quadrature
    /*!
     *  Function that is called to perform the numerical quadrature. Sets the result in the quadratureResult_ local
     *  variable.
     */
    void performQuadrature( )
    {
        quadratureResult_ = performTrapezoidalQuadrature( independentVariables_ , dependentVariables_ );
    }

private:

    //! Independent variables.
    std::vector< IndependentVariableType > independentVariables_;

    //! Dependent variables.
    std::vector< DependentVariableType > dependentVariables_;

    //! Computed value of the quadrature, as computed by last call to performQuadrature.
    DependentVariableType quadratureResult_;

};

//! Typede for trapezoidal quadrature with double (in)dependent variables.
typedef std::shared_ptr< TrapezoidNumericalQuadrature< double, double > > TrapezoidNumericalIntegratorPointerd;

} // namespace numerical_quadrature

} // namespace tudat

#endif // TUDAT_TRAPEZOIDAL_INTEGRATOR_H

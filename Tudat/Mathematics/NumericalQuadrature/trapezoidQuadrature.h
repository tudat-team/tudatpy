/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TRAPEZOIDAL_INTEGRATOR_H
#define TUDAT_TRAPEZOIDAL_INTEGRATOR_H

#include <vector>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/NumericalQuadrature/numericalQuadrature.h"
namespace tudat
{

namespace numerical_quadrature
{


//! Function to perform numerical quadrature using the trapezoidal method.
/*!
 * Function to perform numerical quadrature using the trapezoidal method.
 * \param independentVariables Values of independent variables at which dependentVariables are given
 * \param dependentVariables Values of function for which the numerical quadrature is to be computed, given at
 * independentVariables.
 * \return Numerical quadrature (integral) of the data provided as input
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
     * Constructor
     * \param independentVariables Values of independent variables at which dependentVariables are given
     * \param dependentVariables Values of function for which the numerical quadrature is to be computed, given at
     * independentVariables.
     */
    TrapezoidNumericalQuadrature( const std::vector< IndependentVariableType >& independentVariables,
                                  const std::vector< DependentVariableType >& dependentVariables)
    {
        independentVariables_ = independentVariables;
        dependentVariables_ = dependentVariables;
        performQuadrature( );
    }

    //! Function to reset the (in)dependent variable values.
    /*!
     * Function to reset the (in)dependent variable values.
     * \param independentVariables Values of independent variables at which dependentVariables are given
     * \param dependentVariables Values of function for which the numerical quadrature is to be computed, given at
     * independentVariables.
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
     * Function that is called to perform the numerical quadrature. Sets the result in the quadratureResult_ local
     * variable.
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
typedef boost::shared_ptr< TrapezoidNumericalQuadrature< double, double > > TrapezoidNumericalIntegratorPointerd;

} // namespace numerical_quadrature

} // namespace tudat

#endif // TUDAT_TRAPEZOIDAL_INTEGRATOR_H

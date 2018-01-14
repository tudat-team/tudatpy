/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2002.
 *      Fornberg, B., "Generation of Finite Difference Formulas on Arbitrarily Spaced Grids",
 *          Mathematics of Computation, October 1988.
 *
 */

#ifndef TUDAT_NUMERICAL_DERIVATIVE_H
#define TUDAT_NUMERICAL_DERIVATIVE_H

#include <map>

#include <boost/function.hpp>

#include <Eigen/Core>

namespace tudat
{

namespace numerical_derivatives
{

//! Enum with list of available orders.
/*!
 * List of available orders: 2, 4, 6, 8.
 */
enum CentralDifferenceOrders
{
    order2 = 2,
    order4 = 4,
    order6 = 6,
    order8 = 8
};

//! Get coefficients of a certain order for central difference numerical derivatives.
/*!
 * Get coefficients of a certain order for central difference numerical derivatives as a map
 * of position versus weight. Coefficients are from (Fornberg, B., 1988).
 * \param order The order of the coeffients.
 * \return Map of position versus weight.
 */
const std::map< int, double >& getCentralDifferenceCoefficients( CentralDifferenceOrders order );

//! Compute a numerical derivative using a central difference method.
/*!
 * Computes a specified numerical derivative using a central difference method for a vector
 * function with vector output. The implemented orders are 2nd, 4th and 8th.
 * This function needs to be fully implemented in the header file because it is a template
 * function.
 *
 * \param input Input vector.
 * \param derivativeIndex The index of the entry of the input vector with respect to which to take
 *          the derivative.
 * \param function The function to call with signature
 *          void( const VectorXd& input, VectorXd& result ).
 * \param minimumStep The absolute minimum step size to take. By default 2^-13.
 * \param relativeStepSize The relative step size to take. By default 2^-26.
 * \param order The order of the algorithm to use.
 * \return Numerical derivative calculated from input
 */
template < typename InputType, typename ResultType >
ResultType computeCentralDifference( const InputType& input, const int derivativeIndex,
                                     const boost::function< ResultType( const InputType& ) >& function,
                                     double minimumStep = 0.0, double relativeStepSize = 0.0,
                                     CentralDifferenceOrders order = order2 )
{
    const std::map< int, double >& coefficients = getCentralDifferenceCoefficients( order );

    if ( minimumStep == 0.0 )
    {
        // Set the minimum step to a fourth of the amount of significant
        // digits in a double precision floating point.
        minimumStep = std::pow( 2.0, -13 );
    }

    if ( relativeStepSize == 0.0 )
    {
        // Set the relative step to half of the amount of significant digits
        // in a double precision floating point.
        relativeStepSize = std::pow( 2.0, -26 );
    }

    // Ensure proper rounding by storing the step in a temporary volatile, see
    // (Press W.H., et al., 2002).
    const volatile double temporaryVariable = input( derivativeIndex ) +
            std::max( minimumStep, std::abs( relativeStepSize * input( derivativeIndex ) ) );
    const double realStepSize = temporaryVariable - input( derivativeIndex );

    ResultType result;

    // Compute the numerical derivative.
    for ( std::map< int, double >::const_iterator coefficient = coefficients.begin( );
          coefficient != coefficients.end( ); coefficient++ )
    {
        // Generate deviated input.
        InputType deviatedInput( input );
        deviatedInput( derivativeIndex ) += coefficient->first * realStepSize;

        // Evaluate the function.
        ResultType deviatedResult = function( deviatedInput );

        // Store the result.
        if ( result.size( ) == 0 )
        {
            // Initialize the result.
            result = ResultType::Zero( deviatedResult.rows( ), deviatedResult.cols( ) );
        }

        result += deviatedResult * ( coefficient->second / realStepSize );
    }

    return result;
}

//! Compute a full numerical derivative using a central difference method.
/*!
 * Computes a Jacobian numerically using a central difference method for a vector
 * function with vector output. The implemented orders are 2nd, 4th and 8th.
 *
 * \param input Input vector
 * \param function The function to call with signature
 *          void( const VectorXd& input, VectorXd& result ).
 * \param minimumStep The absolute minimum step size to take. By default 2^-13.
 * \param relativeStepSize The relative step size to take. By default 2^-26.
 * \param order The order of the algorithm to use. Will yield an assertion failure if not 2 or 4.
 * \return Numerical derivative calculated from input
 */
Eigen::MatrixXd computeCentralDifference( const Eigen::VectorXd& input, const boost::function<
                                          Eigen::VectorXd( const Eigen::VectorXd& ) >& function,
                                          double minimumStep = 0.0, double relativeStepSize = 0.0,
                                          CentralDifferenceOrders order = order2 );

template< typename DependentVariableType, typename IndependentVariableType >
DependentVariableType computeCentralDifference(
        const boost::function< DependentVariableType( const IndependentVariableType ) >& dependentVariableFunction,
        const IndependentVariableType nominalIndependentVariable,
        const IndependentVariableType independentVariableStepSize,
        CentralDifferenceOrders order = order2 )
{
    const std::map< int, double >& coefficients = getCentralDifferenceCoefficients( order );

    IndependentVariableType perturbedInput;
    DependentVariableType perturbedOutput;
    DependentVariableType numericalDerivative;

    // Compute the numerical derivative.
    for ( std::map< int, double >::const_iterator coefficientIterator = coefficients.begin( );
          coefficientIterator != coefficients.end( ); coefficientIterator++ )
    {
        // Generate deviated input.
        perturbedInput = nominalIndependentVariable + coefficientIterator->first * independentVariableStepSize;
        perturbedOutput = dependentVariableFunction( perturbedInput );

        if( coefficientIterator == coefficients.begin( ) )
        {
             numericalDerivative = perturbedOutput * ( coefficientIterator->second / independentVariableStepSize );
        }
        else
        {
            numericalDerivative += perturbedOutput * ( coefficientIterator->second / independentVariableStepSize );
        }
    }

    return numericalDerivative;

}

} // namespace numerical_derivatives

} // namespace tudat

#endif // TUDAT_NUMERICAL_DERIVATIVE_H

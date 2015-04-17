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
 *      111214    B. Tong Minh      File added.
 *      120324    K. Kumar          Made minor layout corrections; updated file header to new
 *                                  standard.
 *      120522    E. Heeren         Modified namespace; minor corrections in comments.
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2002.
 *      Fornberg, B., "Generation of Finite Difference Formulas on Arbitrarily Spaced Grids",
 *          Mathematics of Computation, October 1988.
 *
 *    Notes
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
const std::map< int, double >& getCentralDifferenceCoefficients( CentralDifferenceOrders order )
{
    static std::map< CentralDifferenceOrders, std::map< int, double > > coefficients;

    if ( coefficients.empty( ) )
    {
        // Initialize the coefficients map
        coefficients[ order2 ] = std::map< int, double >( );
        coefficients[ order2 ][ -1] = -1.0 / 2.0;
        coefficients[ order2 ][ 1] = 1.0 / 2.0;

        coefficients[ order4 ] = std::map< int, double >( );
        coefficients[ order4 ][ -2] = 1.0 / 12.0;
        coefficients[ order4 ][ -1] = -2.0 / 3.0;
        coefficients[ order4 ][ 1] = 2.0 / 3.0;
        coefficients[ order4 ][ 2] = -1.0 / 12.0;

        coefficients[ order6 ] = std::map< int, double >( );
        coefficients[ order6 ][ -3 ] = -1.0 / 60.0;
        coefficients[ order6 ][ -2 ] = 3.0 / 20.0;
        coefficients[ order6 ][ -1 ] = -3.0 / 4.0;
        coefficients[ order6 ][ 1 ] = 1.0 / 60.0;
        coefficients[ order6 ][ 2 ] = -3.0 / 20.0;
        coefficients[ order6 ][ 3 ] = 3.0 / 4.0;

        coefficients[ order8 ] = std::map< int, double >( );
        coefficients[ order8 ][ -4 ] = 1.0 / 280.0;
        coefficients[ order8 ][ -3 ] = -4.0 / 105.0;
        coefficients[ order8 ][ -2 ] = 1.0 / 5.0;
        coefficients[ order8 ][ -1 ] = -4.0 / 5.0;
        coefficients[ order8 ][ 1 ] = 4.0 / 5.0;
        coefficients[ order8 ][ 2 ] = -1.0 / 5.0;
        coefficients[ order8 ][ 3 ] = 4.0 / 105.0;
        coefficients[ order8 ][ 4 ] = -1.0 / 280.0;
    }

    return coefficients[ order ];
}

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
 * \param stepSize The relative step size to take. By default 2^-26.
 * \param order The order of the algorithm to use.
 * \return Numerical derivative calculated from input
 */
template < typename InputType, typename ResultType >
ResultType computeCentralDifference( const InputType& input, const int derivativeIndex,
                                     const boost::function< ResultType( const InputType& ) >& function,
                                     double minimumStep = 0.0, double stepSize = 0.0,
                                     CentralDifferenceOrders order = order2 )
{
    const std::map< int, double >& coefficients = getCentralDifferenceCoefficients( order );

    if ( minimumStep == 0.0 )
    {
        // Set the minimum step to a fourth of the amount of significant
        // digits in a double precision floating point.
        minimumStep = std::pow( 2.0, -13 );
    }

    if ( stepSize == 0.0 )
    {
        // Set the relative step to half of the amount of significant digits
        // in a double precision floating point.
        stepSize = std::pow( 2.0, -26 );
    }

    // Ensure proper rounding by storing the step in a temporary volatile, see
    // (Press W.H., et al., 2002).
    const volatile double temporaryVariable = input( derivativeIndex ) +
            std::max( minimumStep, std::abs( stepSize * input( derivativeIndex ) ) );
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
 * \param stepSize The relative step size to take. By default 2^-26.
 * \param order The order of the algorithm to use. Will yield an assertion failure if not 2 or 4.
 * \return Numerical derivative calculated from input
 */
Eigen::MatrixXd computeCentralDifference( const Eigen::VectorXd& input, const boost::function<
                                          Eigen::VectorXd( const Eigen::VectorXd& ) >& function,
                                          double minimumStep = 0.0, double stepSize = 0.0,
                                          CentralDifferenceOrders order = order2 )
{
    Eigen::MatrixXd result;

    for ( int derivative = 0; derivative < input.rows( ); derivative++ )
    {
        Eigen::VectorXd partial = computeCentralDifference( input, derivative, function,
                                                            minimumStep, stepSize, order );
        if ( result.size( ) == 0 )
        {
            result = Eigen::MatrixXd( partial.rows( ), input.rows( ) );
        }
        result.col( derivative ) = partial;
    }
    return result;
}

} // namespace numerical_derivatives

} // namespace tudat

#endif // TUDAT_NUMERICAL_DERIVATIVE_H

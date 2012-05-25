/*    Copyright (c) 2011-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      111214    B. Tong Minh      File added.
 *      120324    K. Kumar          Made minor layout corrections; updated file header to new
 *                                  standard.
 *      120522    E. Heeren         Modified namespace; minor corrections in comments.
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
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
namespace basic_mathematics
{
namespace numerical_derivatives
{

//! Enum with list of available orders.
/*!
 * List of available orders: 2, 4, 6, 8.
 */
enum CentralDifferenceOrders
{
    Order2 = 2,
    Order4 = 4,
    Order6 = 6,
    Order8 = 8
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
        coefficients[ Order2 ] = std::map< int, double >( );
        coefficients[ Order2 ][ -1] = -1.0 / 2.0;
        coefficients[ Order2 ][ 1] = 1.0 / 2.0;

        coefficients[ Order4 ] = std::map< int, double >( );
        coefficients[ Order4 ][ -2] = 1.0 / 12.0;
        coefficients[ Order4 ][ -1] = -2.0 / 3.0;
        coefficients[ Order4 ][ 1] = 2.0 / 3.0;
        coefficients[ Order4 ][ 2] = -1.0 / 12.0;

        coefficients[ Order6 ] = std::map< int, double >( );
        coefficients[ Order6 ][ -3 ] = -1.0 / 60.0;
        coefficients[ Order6 ][ -2 ] = 3.0 / 20.0;
        coefficients[ Order6 ][ -1 ] = -3.0 / 4.0;
        coefficients[ Order6 ][ 1 ] = 1.0 / 60.0;
        coefficients[ Order6 ][ 2 ] = -3.0 / 20.0;
        coefficients[ Order6 ][ 3 ] = 3.0 / 4.0;

        coefficients[ Order8 ] = std::map< int, double >( );
        coefficients[ Order8 ][ -4 ] = 1.0 / 280.0;
        coefficients[ Order8 ][ -3 ] = -4.0 / 105.0;
        coefficients[ Order8 ][ -2 ] = 1.0 / 5.0;
        coefficients[ Order8 ][ -1 ] = -4.0 / 5.0;
        coefficients[ Order8 ][ 1 ] = 4.0 / 5.0;
        coefficients[ Order8 ][ 2 ] = -1.0 / 5.0;
        coefficients[ Order8 ][ 3 ] = 4.0 / 105.0;
        coefficients[ Order8 ][ 4 ] = -1.0 / 280.0;
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
 * \param result Vector to store the result to. Caller should initialize this to a zero vector of
 *          proper length.
 * \param function The function to call with signature
 *          void( const VectorXd& input, VectorXd& result ).
 * \param minimumStep The absolute minimum step size to take. By default 2^-13.
 * \param stepSize The relative step size to take. By default 2^-26.
 * \param order The order of the algorithm to use.
 */
template < typename InputType, typename ResultType >
ResultType computeCentralDifference( const InputType& input, const int derivativeIndex,
                                     const boost::function< ResultType( const InputType& ) >& function,
                                     double minimumStep = 0.0, double stepSize = 0.0,
                                     CentralDifferenceOrders order = Order2 )
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
 * \param result Matrix to store the result to. Caller should initialize this to a zero matrix of
 *          proper length.
 * \param function The function to call with signature
 *          void( const VectorXd& input, VectorXd& result ).
 * \param minimumStep The absolute minimum step size to take. By default 2^-13.
 * \param stepSize The relative step size to take. By default 2^-26.
 * \param order The order of the algorithm to use. Will yield an assertion failure if not 2 or 4.
 */
Eigen::MatrixXd computeCentralDifference( const Eigen::VectorXd& input, const boost::function<
                                          Eigen::VectorXd( const Eigen::VectorXd& ) >& function,
                                          double minimumStep = 0.0, double stepSize = 0.0,
                                          CentralDifferenceOrders order = Order2 )
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
} // namespace basic_mathematics
} // namespace tudat

#endif // TUDAT_NUMERICAL_DERIVATIVE_H

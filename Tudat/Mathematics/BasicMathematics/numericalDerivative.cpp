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

#include "Tudat/Mathematics/BasicMathematics/numericalDerivative.h"

namespace tudat
{

namespace numerical_derivatives
{

//! Get coefficients of a certain order for central difference numerical derivatives.
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

Eigen::MatrixXd computeCentralDifference( const Eigen::VectorXd& input, const boost::function<
                                          Eigen::VectorXd( const Eigen::VectorXd& ) >& function,
                                          double minimumStep, double relativeStepSize,
                                          CentralDifferenceOrders order )
{
    Eigen::MatrixXd result;

    for ( int derivative = 0; derivative < input.rows( ); derivative++ )
    {
        Eigen::VectorXd partial = computeCentralDifference( input, derivative, function,
                                                            minimumStep, relativeStepSize, order );
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



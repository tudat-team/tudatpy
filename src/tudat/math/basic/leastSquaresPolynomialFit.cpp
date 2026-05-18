#include "tudat/math/basic/leastSquaresPolynomialFit.h"
#include "tudat/basics/utilities.h"
#include "tudat/math/basic/leastSquaresEstimation.h"

namespace tudat {
namespace linear_algebra {

Eigen::VectorXd evaluatePolynomial( const Eigen::VectorXd& independentValues,
                                    const Eigen::VectorXd& polynomialCoefficients,
                                    const std::vector< double >& polynomialPowers )
{
    Eigen::VectorXd polynomial = Eigen::VectorXd::Zero( independentValues.rows( ) );
    for( unsigned int i = 0; i < polynomialPowers.size( ); i++ )
    {
        polynomial += polynomialCoefficients( i ) * independentValues.array( ).pow( polynomialPowers.at( i ) ).matrix( );
    }
    return polynomial;
}

//! Function to fit a univariate polynomial through a set of data
Eigen::VectorXd getLeastSquaresPolynomialFit( const Eigen::VectorXd& independentValues,
                                              const Eigen::VectorXd& dependentValues,
                                              const std::vector< double >& polynomialPowers )
{
    if( independentValues.rows( ) != dependentValues.rows( ) )
    {
        throw std::runtime_error(
                "Error when doing least squares polynomial fit, size of dependent and independent "
                "variable vectors is not equal." );
    }

    Eigen::MatrixXd designMatrix = Eigen::MatrixXd::Zero( dependentValues.rows( ), polynomialPowers.size( ) );

    // Compute information matrix
    for( int i = 0; i < independentValues.rows( ); i++ )
    {
        for( unsigned int j = 0; j < polynomialPowers.size( ); j++ )
        {
            designMatrix( i, j ) = std::pow( independentValues( i ), polynomialPowers.at( j ) );
        }
    }

    return performLeastSquaresAdjustmentFromDesignMatrix( designMatrix, dependentValues ).first;
}

//! Function to fit a univariate polynomial through a set of data
std::vector< double > getLeastSquaresPolynomialFit( const std::map< double, double >& independentDependentValueMap,
                                                    const std::vector< double >& polynomialPowers )
{
    return utilities::convertEigenVectorToStlVector( getLeastSquaresPolynomialFit(
            utilities::convertStlVectorToEigenVector( utilities::createVectorFromMapKeys( independentDependentValueMap ) ),
            utilities::convertStlVectorToEigenVector( utilities::createVectorFromMapValues( independentDependentValueMap ) ),
            polynomialPowers ) );
}

}}

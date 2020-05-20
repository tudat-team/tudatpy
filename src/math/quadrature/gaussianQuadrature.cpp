#include "tudat/math/quadrature/gaussianQuadrature.h"

namespace tudat
{

namespace numerical_quadrature
{

//! Function to create Gauss quadrature node/weight container
template< typename IndependentVariableType >
std::shared_ptr< GaussQuadratureNodesAndWeights< IndependentVariableType > >
getGaussQuadratureNodesAndWeights( )
{
    return std::make_shared< GaussQuadratureNodesAndWeights< IndependentVariableType > >( );
}

//! Function to create Gauss quadrature node/weight container with long double precision.
template< >
std::shared_ptr< GaussQuadratureNodesAndWeights< long double > >
getGaussQuadratureNodesAndWeights( )
{
    return longDoubleGaussQuadratureNodesAndWeights;
}

//! Function to create Gauss quadrature node/weight container with double precision.
template< >
std::shared_ptr< GaussQuadratureNodesAndWeights< double > >
getGaussQuadratureNodesAndWeights( )
{
    return doubleGaussQuadratureNodesAndWeights;
}

//! Function to create Gauss quadrature node/weight container with float precision.
template< >
std::shared_ptr< GaussQuadratureNodesAndWeights< float > >
getGaussQuadratureNodesAndWeights( )
{
    return floatGaussQuadratureNodesAndWeights;
}


} // namespace numerical_quadrature

} // namespace tudat


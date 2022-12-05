#include "tudat/astro/orbit_determination/podInputOutputTypes.h"

namespace tudat
{

namespace simulation_setup
{

void scaleDesignMatrixWithWeights(
        Eigen::MatrixXd& designMatrix,
        const Eigen::VectorXd& weightsDiagonal )
{
    int numberOfParameters = designMatrix.cols( );
    for( int i = 0; i < weightsDiagonal.rows( ); i++ )
    {
        designMatrix.block( i, 0, 1, numberOfParameters ) *= std::sqrt( weightsDiagonal( i ) );
    }
}


template class EstimationInput< double, double >;
template struct EstimationOutput< double >;

}

}


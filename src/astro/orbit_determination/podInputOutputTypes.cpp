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

Eigen::MatrixXd normaliseUnnormaliseCovarianceMatrix(
        const Eigen::MatrixXd& covarianceMatrix,
        const Eigen::VectorXd& normalisationFactors,
        const bool normalise )
{
    Eigen::MatrixXd modifiedCovarianceMatrix = covarianceMatrix;
    for( int i = 0; i < normalisationFactors.rows( ); i++ )
    {
        for( int j = 0; j < normalisationFactors.rows( ); j++ )
        {
            if ( normalise ) // normalise
            {
                modifiedCovarianceMatrix(i, j) *= normalisationFactors(i) * normalisationFactors(j);
            }
            else // unnormalise
            {
                modifiedCovarianceMatrix( i, j ) /= normalisationFactors( i ) * normalisationFactors( j );
            }
        }
    }
    return modifiedCovarianceMatrix;
}

Eigen::MatrixXd normaliseUnnormaliseInverseCovarianceMatrix(
        Eigen::MatrixXd& inverseCovarianceMatrix,
        Eigen::VectorXd& normalisationFactors,
        const bool normalise )
{
    Eigen::MatrixXd modifiedInverseCovarianceMatrix = inverseCovarianceMatrix;
    for( int i = 0; i < normalisationFactors.rows( ); i++ )
    {
        for( int j = 0; j < normalisationFactors.rows( ); j++ )
        {
            if ( normalise ) // normalise
            {
                modifiedInverseCovarianceMatrix(i, j) /= normalisationFactors(i) * normalisationFactors(j);
            }
            else // unnormalise
            {
                modifiedInverseCovarianceMatrix( i, j ) *= normalisationFactors( i ) * normalisationFactors( j );
            }
        }
    }
    return modifiedInverseCovarianceMatrix;
}

template class CovarianceAnalysisInput< double, double >;
template struct CovarianceAnalysisOutput< double >;

template class EstimationInput< double, double >;
template struct EstimationOutput< double >;

}

}


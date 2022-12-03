#include "tudat/astro/orbit_determination/podInputOutputTypes.h"

namespace tudat
{

namespace simulation_setup
{

void scaleInformationMatrixWithWeights(
        Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& weightsDiagonal )
{
    int numberOfParameters = informationMatrix.cols( );
    for( int i = 0; i < weightsDiagonal.rows( ); i++ )
    {
        informationMatrix.block( i, 0, 1, numberOfParameters ) *= std::sqrt( weightsDiagonal( i ) );
    }
}


template class PodInput< double, double >;
template struct PodOutput< double >;


}

}


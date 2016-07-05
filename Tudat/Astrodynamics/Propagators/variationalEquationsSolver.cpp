#include "Tudat/Astrodynamics/Propagators/variationalEquationsSolver.h"

namespace tudat
{

namespace propagators
{

//! Function to create interpolators for state transition and sensitivity matrices from numerical results.
void createStateTransitionAndSensitivityMatrixInterpolator(
        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >& stateTransitionMatrixInterpolator,
        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >& sensitivityMatrixInterpolator,
        std::vector< std::map< double, Eigen::MatrixXd > >& variationalEquationsSolution,
        const bool clearRawSolution )
{
    // Create interpolator for state transition matrix.
    stateTransitionMatrixInterpolator=
            boost::make_shared< interpolators::LagrangeInterpolator< double, Eigen::MatrixXd > >(
                utilities::createVectorFromMapKeys< Eigen::MatrixXd, double >( variationalEquationsSolution[ 0 ] ),
                utilities::createVectorFromMapValues< Eigen::MatrixXd, double >( variationalEquationsSolution[ 0 ] ), 4 );
    if( clearRawSolution )
    {
        variationalEquationsSolution[ 0 ].clear( );
    }

    // Create interpolator for sensitivity matrix.
    sensitivityMatrixInterpolator =
            boost::make_shared< interpolators::LagrangeInterpolator< double, Eigen::MatrixXd > >(
                utilities::createVectorFromMapKeys< Eigen::MatrixXd, double >( variationalEquationsSolution[ 1 ] ),
                utilities::createVectorFromMapValues< Eigen::MatrixXd, double >( variationalEquationsSolution[ 1 ] ), 4 );
    if( clearRawSolution )
    {
        variationalEquationsSolution[ 1 ].clear( );
    }
}

}

}

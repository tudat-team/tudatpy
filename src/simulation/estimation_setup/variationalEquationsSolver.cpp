/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/estimation_setup/variationalEquationsSolver.h"

namespace tudat
{

namespace propagators
{

template class SingleArcVariationalEquationsSolver< double, double >;
template class MultiArcVariationalEquationsSolver< double, double >;
template class HybridArcVariationalEquationsSolver< double, double >;

//template class VariationalEquationsSolver< double, double >;
//template class VariationalEquationsSolver< long double, double >;
//template class VariationalEquationsSolver< double, Time >;
//template class VariationalEquationsSolver< long double, Time >;

//template class SingleArcVariationalEquationsSolver< double, double >;
//template class SingleArcVariationalEquationsSolver< long double, double >;
//template class SingleArcVariationalEquationsSolver< double, Time >;
//template class SingleArcVariationalEquationsSolver< long double, Time >;

//template class MultiArcVariationalEquationsSolver< double, double >;
//template class MultiArcVariationalEquationsSolver< long double, double >;
//template class MultiArcVariationalEquationsSolver< double, Time >;
//template class MultiArcVariationalEquationsSolver< long double, Time >;

//! Function to create interpolators for state transition and sensitivity matrices from numerical results.
void createStateTransitionAndSensitivityMatrixInterpolator(
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >& stateTransitionMatrixInterpolator,
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >& sensitivityMatrixInterpolator,
        std::map< double, Eigen::MatrixXd >& stateTransitionSolution,
        std::map< double, Eigen::MatrixXd >& sensitivitySolution,
        const bool clearRawSolution )
{
    // Create interpolator for state transition matrix.
    stateTransitionMatrixInterpolator=
            std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::MatrixXd > >(
                utilities::createVectorFromMapKeys< Eigen::MatrixXd, double >( stateTransitionSolution ),
                utilities::createVectorFromMapValues< Eigen::MatrixXd, double >( stateTransitionSolution ), 4 );
//    if( clearRawSolution )
//    {
//        variationalEquationsSolution[ 0 ].clear( );
//    }

    // Create interpolator for sensitivity matrix.
    sensitivityMatrixInterpolator =
            std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::MatrixXd > >(
                utilities::createVectorFromMapKeys< Eigen::MatrixXd, double >( sensitivitySolution ),
                utilities::createVectorFromMapValues< Eigen::MatrixXd, double >( sensitivitySolution ), 4 );

    if( clearRawSolution )
    {
        stateTransitionSolution.clear( );
        sensitivitySolution.clear( );
    }

//    if( clearRawSolution )
//    {
//        variationalEquationsSolution[ 1 ].clear( );
//    }

}



}

}

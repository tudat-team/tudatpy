/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/PropagationSetup/variationalEquationsSolver.h"

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

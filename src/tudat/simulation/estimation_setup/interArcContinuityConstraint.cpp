/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/estimation_setup/interArcContinuityConstraint.h"

#include <Eigen/LU>

namespace tudat
{
namespace simulation_setup
{

int computeConstraintDimension( const Eigen::MatrixXd& constraintWeightMatrix )
{
    // Lari et al. (2021), Eq. (28), defines m_d as the number of scalar jump constraints. Matrix rank is the
    // corresponding count when a general PSD matrix selects arbitrary state-component combinations. Rely on
    // Eigen's scale-aware default rank threshold instead of defining an application-specific eigenvalue cutoff.
    // Validation permits round-off-level asymmetry, so use the same symmetric part as the objective assembly.
    const Eigen::MatrixXd symmetricConstraintWeightMatrix = 0.5 * ( constraintWeightMatrix + constraintWeightMatrix.transpose( ) );
    return symmetricConstraintWeightMatrix.fullPivLu( ).rank( );
}

int computeTotalConstraintDimension( const std::vector< std::shared_ptr< InterArcStateContinuityConstraintSettings > >& constraintSettings )
{
    int totalConstraintDimension = 0;
    for( const auto& settings : constraintSettings )
    {
        if( settings == nullptr )
        {
            throw std::runtime_error( "Error in assembleInterArcContinuityContribution: constraint settings contain a null entry." );
        }
        for( const auto& body : settings->bodies( ) )
        {
            for( std::size_t pairIndex = 0; pairIndex < settings->numberOfPairsForBody( body ); ++pairIndex )
            {
                totalConstraintDimension += computeConstraintDimension( settings->weightMatrixForBodyAndPair( body, pairIndex ) );
            }
        }
    }
    return totalConstraintDimension;
}

}  // namespace simulation_setup
}  // namespace tudat

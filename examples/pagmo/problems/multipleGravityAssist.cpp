/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "multipleGravityAssist.h"

#include "tudat/basics/utilities.h"


//! Descriptive name of the problem
std::string MultipleGravityAssist::get_name() const {
    return "MGA transfer trajectory";
}

//! Get bounds
std::pair<std::vector<double>, std::vector<double> > MultipleGravityAssist::get_bounds() const {

    return { problemBounds_[0], problemBounds_[1] };
}

void MultipleGravityAssist::getDecomposedDecisionVector(
        const Eigen::VectorXd rawDecisionVariables,
        std::vector< double >& currentNodeTimes,
        std::vector< Eigen::VectorXd >& currentLegFreeParameters,
        std::vector< Eigen::VectorXd >& currentNodeFreeParameters ) const
{
    double JD = tudat::physical_constants::JULIAN_DAY;
    currentNodeTimes.push_back( rawDecisionVariables( 0 ) * JD );
    for( unsigned int i = 1; i < numberOfNodes_; i++ )
    {
        currentNodeTimes.push_back(
                    currentNodeTimes.at( i - 1 ) + rawDecisionVariables( i ) * JD );
    }

    for( unsigned int i = 0; i < legParameterIndices_.size( ); i++ )
    {
        currentLegFreeParameters.push_back(
                    rawDecisionVariables.segment( legParameterIndices_.at( i ).first,
                                                  legParameterIndices_.at( i ).second ) );
    }

    for( unsigned int i = 0; i < nodeParameterIndices_.size( ); i++ )
    {
        currentNodeFreeParameters.push_back(
                    rawDecisionVariables.segment( nodeParameterIndices_.at( i ).first,
                                                  nodeParameterIndices_.at( i ).second ) );
    }
}

//! Implementation of the fitness function (return delta-v)
std::vector<double> MultipleGravityAssist::fitness( const std::vector<double> &xv ) const
{
    std::vector< double > currentNodeTimes_;
    std::vector< Eigen::VectorXd > currentLegFreeParameters_;
    std::vector< Eigen::VectorXd > currentNodeFreeParameters_;

    getDecomposedDecisionVector(
                tudat::utilities::convertStlVectorToEigenVector( xv ),
                currentNodeTimes_, currentLegFreeParameters_, currentNodeFreeParameters_ );
    if( transferTrajectory_ != nullptr )
    {
        transferTrajectory_->evaluateTrajectory(
                    currentNodeTimes_, currentLegFreeParameters_,currentNodeFreeParameters_ );
    }
    else
    {
        transferTrajectory_ = createTransferTrajectory(
                    bodyMap_, legSettings_, nodeSettings_, nodeIds_, centralBody_,
                    currentNodeTimes_, currentLegFreeParameters_, currentNodeFreeParameters_);
    }

    return { transferTrajectory_->getTotalDeltaV( ) };
}



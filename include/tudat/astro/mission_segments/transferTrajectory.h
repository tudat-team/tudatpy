/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      Note that the exact implementation of Newton-Raphson as root finder should be updated if
 *      someone would want to use a different root-finding technique.
 *
 *      By default the eccentricity is used as the iteration procedure. This is because in
 *      optimizing a Cassini-like trajectory, the pericenter radius had about 2-4 NaN values in
 *      100000 times the gravity assist calculation. The eccentricity iteration had no NaN values
 *      for a similar run in which 100000 gravity assist calculations were done. Also the
 *      eccentricity seemed to require slightly less iterations (does not necessarily mean it is
 *      faster or more accurate).
 *
 */

#ifndef TUDAT_TRANSFER_TRAJECTORY_H
#define TUDAT_TRANSFER_TRAJECTORY_H

#include <boost/make_shared.hpp>

#include <Eigen/Core>

#include "tudat/astro/mission_segments/transferLeg.h"
#include "tudat/astro/mission_segments/transferNode.h"

namespace tudat
{

namespace mission_segments
{

class TransferTrajectory
{
public:
    TransferTrajectory(
            const std::vector< std::shared_ptr< TransferLeg > > legs,
            const std::vector< std::shared_ptr< TransferNode > > nodes ):
        legs_( legs ), nodes_( nodes ){ }

    void evaluateTrajectory(
            const std::vector< double >& nodeTimes,
            const std::vector< Eigen::VectorXd >& legFreeParameters,
            const std::vector< Eigen::VectorXd >& nodeFreeParameters )
    {
        totalDeltaV_ = 0.0;

        Eigen::VectorXd legTotalParameters;
        for( int i = 0; i < legs_.size( ); i++ )
        {
            getLegTotalParameters(
                        nodeTimes, legFreeParameters.at( i ), i, legTotalParameters );
            legs_.at( i )->updateLegParameters( legTotalParameters );
            totalDeltaV_ += legs_.at( i )->getLegDeltaV( );
        }

        Eigen::VectorXd nodeTotalParameters;
        for( int i = 0; i < nodes_.size( ); i++ )
        {
            getNodeTotalParameters(
                        nodeTimes, nodeFreeParameters.at( i ), i, nodeTotalParameters );
            nodes_.at( i )->updateNodeParameters( nodeTotalParameters );
            totalDeltaV_ += nodes_.at( i )->getNodeDeltaV( );
        }
    }

    double getTotalDeltaV( )
    {
        return totalDeltaV_;
    }

private:

    void getLegTotalParameters(
            const std::vector< double >& nodeTimes,
            const Eigen::VectorXd& legFreeParameters,
            const int legIndex,
            Eigen::VectorXd& legTotalParameters )
    {
        if( legs_.at( legIndex )->getTransferLegType( ) == unpowered_unperturbed_leg )
        {
            legTotalParameters.resize( 2, 1 );
            legTotalParameters( 0 ) = nodeTimes.at( legIndex );
            legTotalParameters( 1 ) = nodeTimes.at( legIndex + 1 );
        }
        else if( legs_.at( legIndex )->getTransferLegType( ) == dsm_position_based_leg )
        {
            legTotalParameters.resize( 6, 1 );
            legTotalParameters( 0 ) = nodeTimes.at( legIndex );
            legTotalParameters( 1 ) = nodeTimes.at( legIndex + 1 );
            legTotalParameters.segment( 2, 4 ) = legFreeParameters;
        }
        else if( legs_.at( legIndex )->getTransferLegType( ) == dsm_velocity_based_leg )
        {
            legTotalParameters.resize( 5, 1 );
            legTotalParameters( 0 ) = nodeTimes.at( legIndex );
            legTotalParameters( 1 ) = nodeTimes.at( legIndex + 1 );
            legTotalParameters.segment( 2, 3 ) = legFreeParameters;
        }
    }

    void getNodeTotalParameters(
            const std::vector< double >& nodeTimes,
            const Eigen::VectorXd& nodeFreeParameters,
            const int nodeIndex,
            Eigen::VectorXd& nodeTotalParameters )
    {
        if( nodes_.at( nodeIndex )->getTransferNodeType( ) == powered_swingby );
        {
            nodeTotalParameters.resize( 1, 1 );
            nodeTotalParameters( 0 ) = nodeTimes.at( i );
        }
        else if( nodes_.at( nodeIndex )->getTransferNodeType( ) == escape_and_departure );
        {
            nodeTotalParameters.resize( 1, 1 );
            nodeTotalParameters( 0 ) = nodeTimes.at( i );
        }
        else if( nodes_.at( nodeIndex )->getTransferNodeType( ) == capture_and_insertion );
        {
            nodeTotalParameters.resize( 1, 1 );
            nodeTotalParameters( 0 ) = nodeTimes.at( i );
        }
    }

    const std::vector< std::shared_ptr< TransferLeg > > legs_;

    const std::vector< std::shared_ptr< TransferNode > > nodes_;

    double totalDeltaV_;
};


} // namespace mission_segments

} // namespace tudat

#endif // TUDAT_TRANSFER_TRAJECTORY_H

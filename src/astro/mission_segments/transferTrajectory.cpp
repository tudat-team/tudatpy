#include "tudat/astro/mission_segments/transferTrajectory.h"

namespace tudat
{

namespace mission_segments
{


void TransferTrajectory::evaluateTrajectory(
        const std::vector< double >& nodeTimes,
        const std::vector< Eigen::VectorXd >& legFreeParameters,
        const std::vector< Eigen::VectorXd >& nodeFreeParameters )
{
    totalDeltaV_ = 0.0;
    totalTimeOfFlight_ = 0.0;

    Eigen::VectorXd legTotalParameters;
    Eigen::VectorXd nodeTotalParameters;

    for( unsigned int i = 0; i < legs_.size( ); i++ )
    {
        if( !nodes_.at( i )->nodeComputesOutgoingVelocity( ) )
        {
            getLegTotalParameters(
                        nodeTimes, legFreeParameters.at( i ), i, legTotalParameters );
            legs_.at( i )->updateLegParameters( legTotalParameters );
            totalDeltaV_ += legs_.at( i )->getLegDeltaV( );
            totalTimeOfFlight_ += legs_.at( i )->getLegTimeOfFlight( );

            getNodeTotalParameters(
                        nodeTimes, nodeFreeParameters.at( i ), i, nodeTotalParameters );
            nodes_.at( i )->updateNodeParameters( nodeTotalParameters );
            totalDeltaV_ += nodes_.at( i )->getNodeDeltaV( );
        }
        else
        {
            getNodeTotalParameters(
                        nodeTimes, nodeFreeParameters.at( i ), i, nodeTotalParameters );
            nodes_.at( i )->updateNodeParameters( nodeTotalParameters );
            totalDeltaV_ += nodes_.at( i )->getNodeDeltaV( );

            getLegTotalParameters(
                        nodeTimes, legFreeParameters.at( i ), i, legTotalParameters );
            legs_.at( i )->updateLegParameters( legTotalParameters );
            totalDeltaV_ += legs_.at( i )->getLegDeltaV( );
            totalTimeOfFlight_ += legs_.at( i )->getLegTimeOfFlight( );
        }
    }

    getNodeTotalParameters(
                nodeTimes, nodeFreeParameters.at( legs_.size( ) ), legs_.size( ), nodeTotalParameters );

    nodes_.at( legs_.size( ) )->updateNodeParameters( nodeTotalParameters );

    totalDeltaV_ += nodes_.at( legs_.size( ) )->getNodeDeltaV( );
    isComputed_ = true;
}

double TransferTrajectory::getTotalDeltaV( )
{
    if( isComputed_ )
    {
        return totalDeltaV_;
    }
    else
    {
        throw std::runtime_error( "Error when getting Delta V for transfer trajectory; transfer parameters not set!" );
    }
}


double TransferTrajectory::getNodeDeltaV( const int nodeIndex )
{
    if( isComputed_ )
    {
        return nodes_.at( nodeIndex )->getNodeDeltaV( );
    }
    else
    {
        throw std::runtime_error( "Error when getting single node Delta V for transfer trajectory; transfer parameters not set!" );
    }
}

double TransferTrajectory::getLegDeltaV( const int legIndex )
{
    if( isComputed_ )
    {
        return legs_.at( legIndex )->getLegDeltaV( );
    }
    else
    {
        throw std::runtime_error( "Error when getting single leg Delta V for transfer trajectory; transfer parameters not set!" );
    }
}

double TransferTrajectory::getTotalTimeOfFlight ( )
{
    if( isComputed_ )
    {
        return totalTimeOfFlight_;
    }
    else
    {
        throw std::runtime_error( "Error when getting Delta V for transfer trajectory; transfer parameters not set!" );
    }
}

std::vector< double > TransferTrajectory::getDeltaVPerNode( )
{
    std::vector< double > deltaVPerNode;
    for( unsigned int i = 0; i < nodes_.size( ); i++ )
    {
        deltaVPerNode.push_back( nodes_.at( i )->getNodeDeltaV( ) );
    }
    return deltaVPerNode;
}

std::vector< double > TransferTrajectory::getDeltaVPerLeg( )
{
    std::vector< double > deltaVPerLeg;
    for( unsigned int i = 0; i < legs_.size( ); i++ )
    {
        deltaVPerLeg.push_back( legs_.at( i )->getLegDeltaV( ) );
    }
    return deltaVPerLeg;
}


void TransferTrajectory::getStatesAlongTrajectoryPerLeg(
        std::vector< std::map< double, Eigen::Vector6d > >& statesAlongTrajectoryPerLeg,
        const int numberOfDataPointsPerLeg )
{
    if( isComputed_ )
    {
        statesAlongTrajectoryPerLeg.clear( );
        statesAlongTrajectoryPerLeg.resize( legs_.size( ) );

        for( unsigned int i = 0; i < legs_.size( ); i++ )
        {
            legs_.at( i )->getStatesAlongTrajectory( statesAlongTrajectoryPerLeg[ i ], numberOfDataPointsPerLeg );
        }
    }
    else
    {
        throw std::runtime_error( "Error when getting states on transfer trajectory; transfer parameters not set!" );
    }
}

void TransferTrajectory::getLegTotalParameters(
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
        legTotalParameters.resize( 3, 1 );
        legTotalParameters( 0 ) = nodeTimes.at( legIndex );
        legTotalParameters( 1 ) = nodeTimes.at( legIndex + 1 );
        legTotalParameters( 2 ) = legFreeParameters( 0 );
    }

}

void TransferTrajectory::getNodeTotalParameters(
        const std::vector< double >& nodeTimes,
        const Eigen::VectorXd& nodeFreeParameters,
        const int nodeIndex,
        Eigen::VectorXd& nodeTotalParameters )
{

    if( nodes_.at( nodeIndex )->getTransferNodeType( ) == swingby )
    {
        if( !nodes_.at( nodeIndex )->nodeComputesOutgoingVelocity( ) )
        {
            nodeTotalParameters.resize( 1, 1 );
            nodeTotalParameters( 0 ) = nodeTimes.at( nodeIndex );
        }
        else
        {
            nodeTotalParameters.resize( 4, 1 );
            nodeTotalParameters( 0 ) = nodeTimes.at( nodeIndex );
            nodeTotalParameters.segment( 1, 3 ) = nodeFreeParameters;
        }
    }
    else if( nodes_.at( nodeIndex )->getTransferNodeType( ) == escape_and_departure )
    {
        if( !nodes_.at( nodeIndex )->nodeComputesOutgoingVelocity( ) )
        {
            nodeTotalParameters.resize( 1, 1 );
            nodeTotalParameters( 0 ) = nodeTimes.at( nodeIndex );
        }
        else
        {
            nodeTotalParameters.resize( 4, 1 );
            nodeTotalParameters( 0 ) = nodeTimes.at( nodeIndex );
            nodeTotalParameters.segment( 1, 3 ) = nodeFreeParameters;
        }
    }
    else if( nodes_.at( nodeIndex )->getTransferNodeType( ) == capture_and_insertion )
    {
        if( !nodes_.at( nodeIndex )->nodeComputesOutgoingVelocity( ) )
        {
            nodeTotalParameters.resize( 1, 1 );
            nodeTotalParameters( 0 ) = nodeTimes.at( nodeIndex );
        }
        else
        {
            throw std::runtime_error( "Error when getting input parameters for capture_and_insertion, node cannot compute output velocity" );
        }
    }
}

} // namespace mission_segments

} // namespace tudat


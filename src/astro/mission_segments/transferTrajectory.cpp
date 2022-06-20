#include "tudat/astro/mission_segments/transferTrajectory.h"
#include "tudat/astro/low_thrust/shape_based/hodographicShapingLeg.h"

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

    std::vector< bool > legEvaluated ( legs_.size( ), false );
    std::vector< bool > nodeEvaluated ( nodes_.size( ), false );

    Eigen::VectorXd legTotalParameters;
    Eigen::VectorXd nodeTotalParameters;

    // Loop over nodes and legs until all are defined
    unsigned int iteration = 0;
    while ( ( std::find(legEvaluated.begin(), legEvaluated.end(), false) != legEvaluated.end() ) ||
        ( std::find(nodeEvaluated.begin(), nodeEvaluated.end(), false) != nodeEvaluated.end() ) )
    {
        ++iteration;

        // First node
        if ( !nodeEvaluated.at( 0 ) && ( nodes_.at( 0 )->nodeComputesOutgoingVelocity( ) || legEvaluated.at( 0 ) ) )
        {
            getNodeTotalParameters( nodeTimes, nodeFreeParameters.at( 0 ), 0, nodeTotalParameters );
            nodes_.at( 0 )->updateNodeParameters( nodeTotalParameters );
            nodeEvaluated.at( 0 ) = true;
            totalDeltaV_ += nodes_.at( 0 )->getNodeDeltaV( );
        }

        // Legs and intermediate nodes
        for( unsigned int i = 0; i < legs_.size( ); i++ )
        {
            // Evaluate leg i
            if ( !legEvaluated.at( i ) && (!nodes_.at( i )->nodeComputesOutgoingVelocity( ) || nodeEvaluated.at( i ) ) &&
                 (!nodes_.at( i+1 )->nodeComputesIncomingVelocity( ) || nodeEvaluated.at( i+1 ) ) )
            {
                getLegTotalParameters( nodeTimes, legFreeParameters.at( i ), i, legTotalParameters );
                legs_.at( i )->updateLegParameters( legTotalParameters );
                legEvaluated.at( i ) = true;
                totalDeltaV_ += legs_.at( i )->getLegDeltaV( );
                totalTimeOfFlight_ += legs_.at( i )->getLegTimeOfFlight( );
            }

            // Evaluate node i+1 (as long as it isn't the last node)
            if ( i < legs_.size( ) - 1 )
            {
                if ( !nodeEvaluated.at( i+1 ) && (nodes_.at( i+1 )->nodeComputesIncomingVelocity( ) || legEvaluated.at(i) ) &&
                     (nodes_.at( i+1 )->nodeComputesOutgoingVelocity( ) || legEvaluated.at(i+1) ) )
                {
                    getNodeTotalParameters( nodeTimes, nodeFreeParameters.at( i+1 ), i+1, nodeTotalParameters );
                    nodes_.at( i+1 )->updateNodeParameters( nodeTotalParameters );
                    nodeEvaluated.at( i+1 ) = true;
                    totalDeltaV_ += nodes_.at( i+1 )->getNodeDeltaV( );
                }
            }
        }

        // Last node
        if ( !nodeEvaluated.at( legs_.size( ) ) && ( nodes_.at( legs_.size( ) )->nodeComputesIncomingVelocity( ) ||
            legEvaluated.at( legs_.size( )-1 ) ) )
        {
            getNodeTotalParameters(nodeTimes, nodeFreeParameters.at( legs_.size( ) ), legs_.size( ), nodeTotalParameters );
            nodes_.at( legs_.size( ) )->updateNodeParameters( nodeTotalParameters );
            nodeEvaluated.at( legs_.size( ) ) = true;
            totalDeltaV_ += nodes_.at( legs_.size( ) )->getNodeDeltaV( );
        }

        if (iteration > legs_.size( ) + nodes_.size() )
        {
            throw std::runtime_error( "evaluateTrajectory used more than maximum possible number of iterations." );
        }

    }

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
        throw std::runtime_error( "Error when getting time of flight for transfer trajectory; transfer parameters not set!" );
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

std::vector< std::map< double, Eigen::Vector6d > > TransferTrajectory::getStatesAlongTrajectoryPerLeg(
        const int numberOfDataPointsPerLeg )
{
    if( isComputed_ )
    {
        std::vector< std::map< double, Eigen::Vector6d > > statesAlongTrajectoryPerLeg;
        getStatesAlongTrajectoryPerLeg( statesAlongTrajectoryPerLeg, numberOfDataPointsPerLeg );
        return statesAlongTrajectoryPerLeg;
    }
    else
    {
        throw std::runtime_error( "Error when getting states on transfer trajectory; transfer parameters not set!" );
    }
}

void TransferTrajectory::getStatesAlongTrajectory( std::map< double, Eigen::Vector6d >& statesAlongTrajectory,
                                                   const int numberOfDataPointsPerLeg )
{
    if( isComputed_ )
    {
        std::vector< std::map< double, Eigen::Vector6d > > statesAlongTrajectoryPerLeg;
        getStatesAlongTrajectoryPerLeg( statesAlongTrajectoryPerLeg, numberOfDataPointsPerLeg );
        statesAlongTrajectory = statesAlongTrajectoryPerLeg.at( 0 );
        for( unsigned int i = 0; i < statesAlongTrajectoryPerLeg.size( ); i++ )
        {
            statesAlongTrajectory.insert( statesAlongTrajectoryPerLeg.at( i ).begin( ),
                                          statesAlongTrajectoryPerLeg.at( i ).end( ) );
        }
    }
    else
    {
        throw std::runtime_error( "Error when getting states on transfer trajectory; transfer parameters not set!" );
    }
}

std::map< double, Eigen::Vector6d > TransferTrajectory::getStatesAlongTrajectory( const int numberOfDataPointsPerLeg )
{
    if( isComputed_ )
    {
        std::map< double, Eigen::Vector6d > statesAlongTrajectory;
        getStatesAlongTrajectory( statesAlongTrajectory, numberOfDataPointsPerLeg );
        return statesAlongTrajectory;
    }
    else
    {
        throw std::runtime_error( "Error when getting states on transfer trajectory; transfer parameters not set!" );
    }
}

void TransferTrajectory::getInertialThrustAccelerationsAlongTrajectoryPerLeg(
        std::vector< std::map< double, Eigen::Vector3d > >& thrustAccelerationsAlongTrajectoryPerLeg,
        const int numberOfDataPointsPerLeg )
{
    if( isComputed_ )
    {
        thrustAccelerationsAlongTrajectoryPerLeg.clear( );
        thrustAccelerationsAlongTrajectoryPerLeg.resize( legs_.size( ) );

        for( unsigned int i = 0; i < legs_.size( ); i++ )
        {
            legs_.at( i )->getThrustAccelerationsAlongTrajectory( thrustAccelerationsAlongTrajectoryPerLeg[ i ],
                                                                  numberOfDataPointsPerLeg );
        }
    }
    else
    {
        throw std::runtime_error( "Error when getting thrust accelerations on transfer trajectory; transfer parameters not set!" );
    }
}

std::vector< std::map< double, Eigen::Vector3d > > TransferTrajectory::getInertialThrustAccelerationAlongTrajectoryPerLeg(
        const int numberOfDataPointsPerLeg )
{
    if( isComputed_ )
    {
        std::vector< std::map< double, Eigen::Vector3d > > thrustAccelerationsAlongTrajectoryPerLeg;
        getInertialThrustAccelerationsAlongTrajectoryPerLeg(thrustAccelerationsAlongTrajectoryPerLeg,
                                                            numberOfDataPointsPerLeg);
        return thrustAccelerationsAlongTrajectoryPerLeg;
    }
    else
    {
        throw std::runtime_error( "Error when getting thrust accelerations on transfer trajectory; transfer parameters not set!" );
    }
}

void TransferTrajectory::getInertialThrustAccelerationsAlongTrajectory(
        std::map< double, Eigen::Vector3d >& thrustAccelerationsAlongTrajectory,
        const int numberOfDataPointsPerLeg )
{
    if( isComputed_ )
    {
        std::vector< std::map< double, Eigen::Vector3d > > thrustAccelerationsAlongTrajectoryPerLeg;
        getInertialThrustAccelerationsAlongTrajectoryPerLeg(thrustAccelerationsAlongTrajectoryPerLeg,
                                                            numberOfDataPointsPerLeg);
        thrustAccelerationsAlongTrajectory = thrustAccelerationsAlongTrajectoryPerLeg.at( 0 );
        for( unsigned int i = 0; i < thrustAccelerationsAlongTrajectoryPerLeg.size( ); i++ )
        {
            thrustAccelerationsAlongTrajectory.insert( thrustAccelerationsAlongTrajectoryPerLeg.at( i ).begin( ),
                                                       thrustAccelerationsAlongTrajectoryPerLeg.at( i ).end( ) );
        }
    }
    else
    {
        throw std::runtime_error( "Error when getting thrust accelerations on transfer trajectory; transfer parameters not set!" );
    }
}

std::map< double, Eigen::Vector3d > TransferTrajectory::getInertialThrustAccelerationsAlongTrajectory(
        const int numberOfDataPointsPerLeg )
{
    if( isComputed_ )
    {
        std::map< double, Eigen::Vector3d > thrustAccelerationsAlongTrajectory;
        getInertialThrustAccelerationsAlongTrajectory(thrustAccelerationsAlongTrajectory, numberOfDataPointsPerLeg);
        return thrustAccelerationsAlongTrajectory;
    }
    else
    {
        throw std::runtime_error( "Error when getting thrust accelerations on transfer trajectory; transfer parameters not set!" );
    }
}

std::map< double, Eigen::Vector3d > TransferTrajectory::getRswThrustAccelerationsAlongTrajectory( const int numberOfDataPointsPerLeg )
{
    if( isComputed_ )
    {
        // Get inertial acceleration
        std::map< double, Eigen::Vector3d > inertialThrustAccelerationsAlongTrajectory;
        getInertialThrustAccelerationsAlongTrajectory(inertialThrustAccelerationsAlongTrajectory,
                                                      numberOfDataPointsPerLeg);

        // Get inertial state
        std::map< double, Eigen::Vector6d > statesAlongTrajectory;
        getStatesAlongTrajectory(statesAlongTrajectory, numberOfDataPointsPerLeg );

        std::map< double, Eigen::Vector3d > rswThrustAccelerationsAlongTrajectory;

        for( std::map< double, Eigen::Vector3d >::iterator it = inertialThrustAccelerationsAlongTrajectory.begin();
            it != inertialThrustAccelerationsAlongTrajectory.end(); ++it )
        {
            double time = it->first;
            rswThrustAccelerationsAlongTrajectory[time] =
                    reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrix(statesAlongTrajectory[time]) *
                    inertialThrustAccelerationsAlongTrajectory[time];
        }

        return rswThrustAccelerationsAlongTrajectory;
    }
    else
    {
        throw std::runtime_error( "Error when getting thrust accelerations on transfer trajectory; transfer parameters not set!" );
    }
}

std::map< double, Eigen::Vector3d > TransferTrajectory::getTnwThrustAccelerationsAlongTrajectory( const int numberOfDataPointsPerLeg )
{
    if( isComputed_ )
    {
        // Get inertial acceleration
        std::map< double, Eigen::Vector3d > inertialThrustAccelerationsAlongTrajectory;
        getInertialThrustAccelerationsAlongTrajectory(inertialThrustAccelerationsAlongTrajectory,
                                                      numberOfDataPointsPerLeg);

        // Get inertial state
        std::map< double, Eigen::Vector6d > statesAlongTrajectory;
        getStatesAlongTrajectory(statesAlongTrajectory, numberOfDataPointsPerLeg );

        std::map< double, Eigen::Vector3d > tnwThrustAccelerationsAlongTrajectory;

        for( std::map< double, Eigen::Vector3d >::iterator it = inertialThrustAccelerationsAlongTrajectory.begin();
            it != inertialThrustAccelerationsAlongTrajectory.end(); ++it )
        {
            double time = it->first;
            tnwThrustAccelerationsAlongTrajectory[time] =
                    reference_frames::getInertialToTnwRotation(statesAlongTrajectory[time]) *
                    inertialThrustAccelerationsAlongTrajectory[time];
        }

        return tnwThrustAccelerationsAlongTrajectory;
    }
    else
    {
        throw std::runtime_error( "Error when getting thrust accelerations on transfer trajectory; transfer parameters not set!" );
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
    else if( legs_.at( legIndex )->getTransferLegType( ) == spherical_shaping_low_thrust_leg )
    {
        legTotalParameters.resize( 3, 1 );
        legTotalParameters( 0 ) = nodeTimes.at( legIndex );
        legTotalParameters( 1 ) = nodeTimes.at( legIndex + 1 );
        legTotalParameters( 2 ) = legFreeParameters( 0 );
    }
    else if ( legs_.at( legIndex )->getTransferLegType( ) == hodographic_low_thrust_leg )
    {
        std::shared_ptr< shape_based_methods::HodographicShapingLeg > hodographicShapingLeg =
                std::dynamic_pointer_cast< shape_based_methods::HodographicShapingLeg >( legs_.at( legIndex ) );
        if (hodographicShapingLeg == nullptr)
        {
            throw std::runtime_error("Error getting leg parameters, hodographic leg shaping type is invalid");
        }
        const int legFreeVelocityShapingParameters = hodographicShapingLeg->getNumberOfFreeCoefficients();

        legTotalParameters.resize( 3 + legFreeVelocityShapingParameters, 1 );
        legTotalParameters( 0 ) = nodeTimes.at( legIndex );
        legTotalParameters( 1 ) = nodeTimes.at( legIndex + 1 );
        legTotalParameters.segment( 2, 1 + legFreeVelocityShapingParameters) = legFreeParameters;
    }
    else
    {
        throw std::runtime_error( "Error when getting leg parameters, leg type not recognized" );
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
        if( !nodes_.at( nodeIndex )->nodeComputesOutgoingVelocity( ) && !nodes_.at( nodeIndex )->nodeComputesIncomingVelocity( ) )
        {
            nodeTotalParameters.resize( 1, 1 );
            nodeTotalParameters( 0 ) = nodeTimes.at( nodeIndex );
        }
        else if ( nodes_.at( nodeIndex )->nodeComputesOutgoingVelocity( ) && !nodes_.at( nodeIndex )->nodeComputesIncomingVelocity( ) )
        {
            nodeTotalParameters.resize( 4, 1 );
            nodeTotalParameters( 0 ) = nodeTimes.at( nodeIndex );
            nodeTotalParameters.segment( 1, 3 ) = nodeFreeParameters;
        }
        else if ( nodes_.at( nodeIndex )->nodeComputesOutgoingVelocity( ) && nodes_.at( nodeIndex )->nodeComputesIncomingVelocity( ) )
        {
            nodeTotalParameters.resize( 7, 1 );
            nodeTotalParameters( 0 ) = nodeTimes.at( nodeIndex );
            nodeTotalParameters.segment( 1, 6 ) = nodeFreeParameters;
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
        if( !nodes_.at( nodeIndex )->nodeComputesIncomingVelocity( ) )
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
    else
    {
        throw std::runtime_error( "Error when getting node parameters, node type not recognized" );
    }
}

} // namespace mission_segments

} // namespace tudat


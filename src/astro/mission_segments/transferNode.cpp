#include "tudat/astro/mission_segments/transferNode.h"

namespace tudat
{
namespace mission_segments
{


TransferNode::TransferNode(
        const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
        const TransferNodeTypes nodeType ):
    nodeEphemeris_( nodeEphemeris ), nodeType_( nodeType ), nodeParameters_( Eigen::VectorXd::Zero( 0 ) ){ }

void TransferNode::updateNodeParameters( const Eigen::VectorXd nodeParameters )
{
    nodeParameters_ = nodeParameters;
    computeNode( );
}

double TransferNode::getNodeDeltaV( )
{
    return totalNodeDeltaV_;
}

TransferNodeTypes TransferNode::getTransferNodeType( )
{
    return nodeType_;
}

Eigen::Vector3d TransferNode::getIncomingVelocity( )
{
    return incomingVelocity_;
}

Eigen::Vector3d TransferNode::getOutgoingVelocity( )
{
    return outgoingVelocity_;
}


void TransferNode::updateNodeState( const double nodeTime )
{
    nodeState_ = nodeEphemeris_->getCartesianState( nodeTime );
}

DepartureWithFixedOutgoingVelocityNode::DepartureWithFixedOutgoingVelocityNode(
        const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
        const double centralBodyGravitationalParameter,
        const double departureSemiMajorAxis,
        const double departureEccentricity,
        const std::function< Eigen::Vector3d( ) > outgoingVelocityFunction ):
    TransferNode( nodeEphemeris, escape_and_departure ),
    centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
    departureSemiMajorAxis_( departureSemiMajorAxis ),
    departureEccentricity_( departureEccentricity ),
    outgoingVelocityFunction_( outgoingVelocityFunction )
{}

Eigen::Vector3d DepartureWithFixedOutgoingVelocityNode::getIncomingVelocity( )
{
    throw std::runtime_error( "Error, no incoming velocity can be given for escape_and_departure leg" );
}

bool DepartureWithFixedOutgoingVelocityNode::nodeComputesOutgoingVelocity( )
{
    return false;
}

void DepartureWithFixedOutgoingVelocityNode::computeNode( )
{
    if( nodeParameters_.rows( ) != 1 )
    {
        throw std::runtime_error( "Error when computing DepartureWithFixedOutgoingVelocityNode, incorrect input size" );
    }

    nodeTime_ = nodeParameters_( 0 );
    updateNodeState( nodeTime_ );

    outgoingVelocity_ = outgoingVelocityFunction_( );

    totalNodeDeltaV_ = mission_segments::computeEscapeOrCaptureDeltaV(
                centralBodyGravitationalParameter_,
                departureSemiMajorAxis_, departureEccentricity_,
                ( outgoingVelocity_ - nodeState_.segment< 3 >( 3 ) ).norm( ) );
}





DepartureWithFreeOutgoingVelocityNode::DepartureWithFreeOutgoingVelocityNode(
        const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
        const double centralBodyGravitationalParameter,
        const double departureSemiMajorAxis,
        const double departureEccentricity ):
    TransferNode( nodeEphemeris, escape_and_departure ),
    centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
    departureSemiMajorAxis_( departureSemiMajorAxis ),
    departureEccentricity_( departureEccentricity )
{ }

Eigen::Vector3d DepartureWithFreeOutgoingVelocityNode::getIncomingVelocity( )
{
    throw std::runtime_error( "Error, no incoming velocity can be given for escape_and_departure leg" );
}

bool DepartureWithFreeOutgoingVelocityNode::nodeComputesOutgoingVelocity( )
{
    return true;
}

void DepartureWithFreeOutgoingVelocityNode::computeNode( )
{
    if( nodeParameters_.rows( ) != 4 )
    {
        throw std::runtime_error( "Error when computing DepartureWithFreeOutgoingVelocityNode, incorrect input size" );
    }

    nodeTime_ = nodeParameters_( 0 );
    excessVelocityMagnitude_ = nodeParameters_( 1 );
    excessVelocityInPlaneAngle_ = nodeParameters_( 2 );
    excessVelocityOutOfPlaneAngle_ = nodeParameters_( 3 );

    updateNodeState( nodeTime_ );

    // Calculate unit vectors as described in [Vinko and Izzo, 2008].
    Eigen::Vector3d nodeVelocity = nodeState_.segment< 3 >( 3 );
    const Eigen::Vector3d unitVector1 = nodeVelocity.normalized( );
    const Eigen::Vector3d unitVector3 = ( nodeState_.segment< 3 >( 0 ).cross( nodeVelocity ) ).normalized( );
    const Eigen::Vector3d unitVector2 = unitVector3.cross( unitVector1 );


    // Calculate the velocity after departure as described in [Vinko and Izzo, 2008].
    outgoingVelocity_ = nodeVelocity +
            excessVelocityMagnitude_ * std::cos( excessVelocityInPlaneAngle_ ) *
            std::cos( excessVelocityOutOfPlaneAngle_ ) * unitVector1 +
            excessVelocityMagnitude_ * std::sin( excessVelocityInPlaneAngle_ ) *
            std::cos( excessVelocityOutOfPlaneAngle_ ) * unitVector2 +
            excessVelocityMagnitude_ * std::sin( excessVelocityOutOfPlaneAngle_ ) * unitVector3;

    totalNodeDeltaV_ = mission_segments::computeEscapeOrCaptureDeltaV(
                centralBodyGravitationalParameter_,
                departureSemiMajorAxis_, departureEccentricity_,
                ( outgoingVelocity_ - nodeVelocity ).norm( ) );

}


CaptureAndInsertionNode::CaptureAndInsertionNode(
        const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
        const double centralBodyGravitationalParameter,
        const double captureSemiMajorAxis,
        const double captureEccentricity,
        const std::function< Eigen::Vector3d( ) > incomingVelocityFunction ):
    TransferNode( nodeEphemeris, capture_and_insertion ),
    centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
    captureSemiMajorAxis_( captureSemiMajorAxis ),
    captureEccentricity_( captureEccentricity ),
    incomingVelocityFunction_( incomingVelocityFunction )
{ }

Eigen::Vector3d CaptureAndInsertionNode::getOutgoingVelocity( )
{
    throw std::runtime_error( "Error, no outgoing velocity can be given for capture_and_insertion leg" );
}

bool CaptureAndInsertionNode::nodeComputesOutgoingVelocity( )
{
    return false;
}

void CaptureAndInsertionNode::computeNode( )
{
    if( nodeParameters_.rows( ) != 1 )
    {
        throw std::runtime_error( "Error when computing CaptureAndInsertionNode, incorrect input size" );
    }

    nodeTime_ = nodeParameters_( 0 );

    updateNodeState( nodeTime_ );

    incomingVelocity_ = incomingVelocityFunction_( );

    totalNodeDeltaV_ = mission_segments::computeEscapeOrCaptureDeltaV(
                centralBodyGravitationalParameter_,
                captureSemiMajorAxis_, captureEccentricity_,
                ( incomingVelocity_ - nodeState_.segment< 3 >( 3 ) ).norm( ) );
}


SwingbyWithFixedOutgoingVelocity::SwingbyWithFixedOutgoingVelocity(
        const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
        const double centralBodyGravitationalParameter,
        const double minimumPeriapsisRadius,
        const std::function< Eigen::Vector3d( ) > incomingVelocityFunction,
        const std::function< Eigen::Vector3d( ) > outgoingVelocityFunction ):
    TransferNode( nodeEphemeris, swingby ),
    centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
    minimumPeriapsisRadius_( minimumPeriapsisRadius ),
    incomingVelocityFunction_( incomingVelocityFunction ),
    outgoingVelocityFunction_( outgoingVelocityFunction )
{ }

bool SwingbyWithFixedOutgoingVelocity::nodeComputesOutgoingVelocity( )
{
    return false;
}


void SwingbyWithFixedOutgoingVelocity::computeNode( )
{
    if( nodeParameters_.rows( ) != 1 )
    {
        throw std::runtime_error( "Error when computing SwingbyWithFixedOutgoingVelocity, incorrect input size" );
    }

    nodeTime_ = nodeParameters_( 0 );

    updateNodeState( nodeTime_ );

    incomingVelocity_ = incomingVelocityFunction_( );
    outgoingVelocity_ = outgoingVelocityFunction_( );

    totalNodeDeltaV_ = calculateGravityAssistDeltaV(
                centralBodyGravitationalParameter_, nodeState_.segment< 3 >( 3 ),
                incomingVelocity_, outgoingVelocity_, minimumPeriapsisRadius_ );

}




SwingbyWithFreeOutgoingVelocity::SwingbyWithFreeOutgoingVelocity(
        const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
        const double centralBodyGravitationalParameter,
        const std::function< Eigen::Vector3d( ) > incomingVelocityFunction ):
    TransferNode( nodeEphemeris, swingby ),
    centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
    incomingVelocityFunction_( incomingVelocityFunction )
{ }

bool SwingbyWithFreeOutgoingVelocity::nodeComputesOutgoingVelocity( )
{
    return true;
}


void SwingbyWithFreeOutgoingVelocity::computeNode( )
{
    if( nodeParameters_.rows( ) != 4 )
    {
        throw std::runtime_error( "Error when computing SwingbyWithFreeOutgoingVelocity, incorrect input size" );
    }

    nodeTime_ = nodeParameters_( 0 );
    periapsisRadius_ = nodeParameters_( 1 );
    outgoingRotationAngle_ = nodeParameters_( 2 );
    swingbyDeltaV_ = nodeParameters_( 3 );

    updateNodeState( nodeTime_ );

    incomingVelocity_ = incomingVelocityFunction_( );

    // Prepare the gravity assist propagator module.
    outgoingVelocity_ = mission_segments::calculatePoweredGravityAssistOutgoingVelocity(
                centralBodyGravitationalParameter_,
                nodeState_.segment< 3 >( 3 ), incomingVelocity_,
                outgoingRotationAngle_, periapsisRadius_, swingbyDeltaV_ );

    totalNodeDeltaV_ = swingbyDeltaV_;

}



} // namespace mission_segments

} // namespace tudat


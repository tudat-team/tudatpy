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

bool DepartureWithFixedOutgoingVelocityNode::nodeComputesIncomingVelocity( )
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

bool DepartureWithFreeOutgoingVelocityNode::nodeComputesIncomingVelocity( )
{
    return false;
}

void DepartureWithFreeOutgoingVelocityNode::computeNode( )
{
    if( nodeParameters_.rows( ) != 4 )
    {
        throw std::runtime_error( "Error when computing DepartureWithFreeOutgoingVelocityNode, incorrect input size" );
    }

    nodeTime_ = nodeParameters_( 0 );
    outgoingExcessVelocityMagnitude_ = nodeParameters_(1 );
    outgoingExcessVelocityInPlaneAngle_ = nodeParameters_(2 );
    outgoingExcessVelocityOutOfPlaneAngle_ = nodeParameters_(3 );

    updateNodeState( nodeTime_ );

    // Calculate unit vectors as described in [Vinko and Izzo, 2008].
    Eigen::Vector3d nodeVelocity = nodeState_.segment< 3 >( 3 );
    const Eigen::Vector3d unitVector1 = nodeVelocity.normalized( );
    const Eigen::Vector3d unitVector3 = ( nodeState_.segment< 3 >( 0 ).cross( nodeVelocity ) ).normalized( );
    const Eigen::Vector3d unitVector2 = unitVector3.cross( unitVector1 );


    // Calculate the outgoing velocity as described in [Vinko and Izzo, 2008].
    outgoingVelocity_ = nodeVelocity +
                        outgoingExcessVelocityMagnitude_ * std::cos(outgoingExcessVelocityInPlaneAngle_ ) *
                        std::cos(outgoingExcessVelocityOutOfPlaneAngle_ ) * unitVector1 +
                        outgoingExcessVelocityMagnitude_ * std::sin(outgoingExcessVelocityInPlaneAngle_ ) *
                        std::cos(outgoingExcessVelocityOutOfPlaneAngle_ ) * unitVector2 +
                        outgoingExcessVelocityMagnitude_ * std::sin(outgoingExcessVelocityOutOfPlaneAngle_ ) * unitVector3;

    totalNodeDeltaV_ = mission_segments::computeEscapeOrCaptureDeltaV(
                centralBodyGravitationalParameter_,
                departureSemiMajorAxis_, departureEccentricity_,
                ( outgoingVelocity_ - nodeVelocity ).norm( ) );
}


CaptureWithFixedIncomingVelocityNode::CaptureWithFixedIncomingVelocityNode(
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

Eigen::Vector3d CaptureWithFixedIncomingVelocityNode::getOutgoingVelocity( )
{
    throw std::runtime_error( "Error, no outgoing velocity can be given for capture_and_insertion node" );
}

bool CaptureWithFixedIncomingVelocityNode::nodeComputesOutgoingVelocity( )
{
    return false;
}

bool CaptureWithFixedIncomingVelocityNode::nodeComputesIncomingVelocity( )
{
    return false;
}

void CaptureWithFixedIncomingVelocityNode::computeNode( )
{
    if( nodeParameters_.rows( ) != 1 )
    {
        throw std::runtime_error( "Error when computing CaptureWithFixedIncomingVelocityNode, incorrect input size" );
    }

    nodeTime_ = nodeParameters_( 0 );

    updateNodeState( nodeTime_ );

    incomingVelocity_ = incomingVelocityFunction_( );

    totalNodeDeltaV_ = mission_segments::computeEscapeOrCaptureDeltaV(
                centralBodyGravitationalParameter_,
                captureSemiMajorAxis_, captureEccentricity_,
                ( incomingVelocity_ - nodeState_.segment< 3 >( 3 ) ).norm( ) );
}


CaptureWithFreeIncomingVelocityNode::CaptureWithFreeIncomingVelocityNode(
        const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
        const double centralBodyGravitationalParameter,
        const double captureSemiMajorAxis,
        const double captureEccentricity):
    TransferNode( nodeEphemeris, capture_and_insertion ),
    centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
    captureSemiMajorAxis_( captureSemiMajorAxis ),
    captureEccentricity_( captureEccentricity )
{ }

Eigen::Vector3d CaptureWithFreeIncomingVelocityNode::getOutgoingVelocity( )
{
    throw std::runtime_error( "Error, no outgoing velocity can be given for capture_and_insertion node" );
}

bool CaptureWithFreeIncomingVelocityNode::nodeComputesOutgoingVelocity( )
{
    return false;
}

bool CaptureWithFreeIncomingVelocityNode::nodeComputesIncomingVelocity( )
{
    return true;
}

void CaptureWithFreeIncomingVelocityNode::computeNode( )
{
    if( nodeParameters_.rows( ) != 4 )
    {
        throw std::runtime_error( "Error when computing CaptureWithFreeIncomingVelocityNode, incorrect input size" );
    }

    nodeTime_ = nodeParameters_( 0 );
    incomingExcessVelocityMagnitude_ = nodeParameters_( 1 );
    incomingExcessVelocityInPlaneAngle_ = nodeParameters_( 2 );
    incomingExcessVelocityOutOfPlaneAngle_ = nodeParameters_( 3 );

    updateNodeState( nodeTime_ );

    // Calculate unit vectors as described in [Vinko and Izzo, 2008].
    Eigen::Vector3d nodeVelocity = nodeState_.segment< 3 >( 3 );
    const Eigen::Vector3d unitVector1 = nodeVelocity.normalized( );
    const Eigen::Vector3d unitVector3 = ( nodeState_.segment< 3 >( 0 ).cross( nodeVelocity ) ).normalized( );
    const Eigen::Vector3d unitVector2 = unitVector3.cross( unitVector1 );

    // Calculate the incoming velocity as described in [Vinko and Izzo, 2008].
    incomingVelocity_ = nodeVelocity +
                        incomingExcessVelocityMagnitude_ * std::cos(incomingExcessVelocityInPlaneAngle_ ) *
                        std::cos(incomingExcessVelocityOutOfPlaneAngle_ ) * unitVector1 +
                        incomingExcessVelocityMagnitude_ * std::sin(incomingExcessVelocityInPlaneAngle_ ) *
                        std::cos(incomingExcessVelocityOutOfPlaneAngle_ ) * unitVector2 +
                        incomingExcessVelocityMagnitude_ * std::sin(incomingExcessVelocityOutOfPlaneAngle_ ) * unitVector3;

    totalNodeDeltaV_ = mission_segments::computeEscapeOrCaptureDeltaV(
                centralBodyGravitationalParameter_,
                captureSemiMajorAxis_, captureEccentricity_,
                ( incomingVelocity_ - nodeState_.segment< 3 >( 3 ) ).norm( ) );
}


SwingbyWithFixedIncomingFixedOutgoingVelocity::SwingbyWithFixedIncomingFixedOutgoingVelocity(
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

bool SwingbyWithFixedIncomingFixedOutgoingVelocity::nodeComputesOutgoingVelocity( )
{
    return false;
}

bool SwingbyWithFixedIncomingFixedOutgoingVelocity::nodeComputesIncomingVelocity( )
{
    return false;
}


void SwingbyWithFixedIncomingFixedOutgoingVelocity::computeNode( )
{
    if( nodeParameters_.rows( ) != 1 )
    {
        throw std::runtime_error( "Error when computing SwingbyWithFixedIncomingFixedOutgoingVelocity, incorrect input size" );
    }

    nodeTime_ = nodeParameters_( 0 );

    updateNodeState( nodeTime_ );

    incomingVelocity_ = incomingVelocityFunction_( );
    outgoingVelocity_ = outgoingVelocityFunction_( );

    totalNodeDeltaV_ = calculateGravityAssistDeltaV(
                centralBodyGravitationalParameter_, nodeState_.segment< 3 >( 3 ),
                incomingVelocity_, outgoingVelocity_, minimumPeriapsisRadius_ );

}


SwingbyWithFixedIncomingFreeOutgoingVelocity::SwingbyWithFixedIncomingFreeOutgoingVelocity(
        const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
        const double centralBodyGravitationalParameter,
        const std::function< Eigen::Vector3d( ) > incomingVelocityFunction ):
    TransferNode( nodeEphemeris, swingby ),
    centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
    incomingVelocityFunction_( incomingVelocityFunction )
{ }

bool SwingbyWithFixedIncomingFreeOutgoingVelocity::nodeComputesOutgoingVelocity( )
{
    return true;
}

bool SwingbyWithFixedIncomingFreeOutgoingVelocity::nodeComputesIncomingVelocity( )
{
    return false;
}

void SwingbyWithFixedIncomingFreeOutgoingVelocity::computeNode( )
{
    if( nodeParameters_.rows( ) != 4 )
    {
        throw std::runtime_error( "Error when computing SwingbyWithFixedIncomingFreeOutgoingVelocity, incorrect input size" );
    }

    nodeTime_ = nodeParameters_( 0 );
    periapsisRadius_ = nodeParameters_( 1 );
    outgoingRotationAngle_ = nodeParameters_( 2 );
    swingbyDeltaV_ = nodeParameters_( 3 );

    updateNodeState( nodeTime_ );

    incomingVelocity_ = incomingVelocityFunction_( );

    // Forward propagate the gravity assist
    outgoingVelocity_ = mission_segments::calculatePoweredGravityAssistOutgoingVelocity(
                centralBodyGravitationalParameter_,
                nodeState_.segment< 3 >( 3 ), incomingVelocity_,
                outgoingRotationAngle_, periapsisRadius_, swingbyDeltaV_ );

    totalNodeDeltaV_ = swingbyDeltaV_;

}


SwingbyWithFreeIncomingFreeOutgoingVelocity::SwingbyWithFreeIncomingFreeOutgoingVelocity(
           const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
           const double centralBodyGravitationalParameter):
       TransferNode( nodeEphemeris, swingby ),
       centralBodyGravitationalParameter_( centralBodyGravitationalParameter )
{ }

bool SwingbyWithFreeIncomingFreeOutgoingVelocity::nodeComputesOutgoingVelocity( )
{
    return true;
}

bool SwingbyWithFreeIncomingFreeOutgoingVelocity::nodeComputesIncomingVelocity( )
{
    return true;
}

void SwingbyWithFreeIncomingFreeOutgoingVelocity::computeNode( )
{
    if( nodeParameters_.rows( ) != 7 )
    {
        throw std::runtime_error( "Error when computing SwingbyWithFreeIncomingFreeOutgoingVelocity, incorrect input size" );
    }

    nodeTime_ = nodeParameters_( 0 );
    incomingExcessVelocityMagnitude_ = nodeParameters_( 1 );
    incomingExcessVelocityInPlaneAngle_ = nodeParameters_( 2 );
    incomingExcessVelocityOutOfPlaneAngle_ = nodeParameters_( 3 );
    periapsisRadius_ = nodeParameters_( 4 );
    outgoingRotationAngle_ = nodeParameters_( 5 );
    swingbyDeltaV_ = nodeParameters_( 6 );

    updateNodeState( nodeTime_ );

    // Calculate unit vectors as described in [Vinko and Izzo, 2008].
    Eigen::Vector3d nodeVelocity = nodeState_.segment< 3 >( 3 );
    const Eigen::Vector3d unitVector1 = nodeVelocity.normalized( );
    const Eigen::Vector3d unitVector3 = ( nodeState_.segment< 3 >( 0 ).cross( nodeVelocity ) ).normalized( );
    const Eigen::Vector3d unitVector2 = unitVector3.cross( unitVector1 );
    // Calculate the incoming velocity as described in [Vinko and Izzo, 2008].
    incomingVelocity_ = nodeVelocity +
                        incomingExcessVelocityMagnitude_ * std::cos(incomingExcessVelocityInPlaneAngle_ ) *
                        std::cos(incomingExcessVelocityOutOfPlaneAngle_ ) * unitVector1 +
                        incomingExcessVelocityMagnitude_ * std::sin(incomingExcessVelocityInPlaneAngle_ ) *
                        std::cos(incomingExcessVelocityOutOfPlaneAngle_ ) * unitVector2 +
                        incomingExcessVelocityMagnitude_ * std::sin(incomingExcessVelocityOutOfPlaneAngle_ ) * unitVector3;

    // Forward propagate the gravity assist
    outgoingVelocity_ = mission_segments::calculatePoweredGravityAssistOutgoingVelocity(
                centralBodyGravitationalParameter_,
                nodeState_.segment< 3 >( 3 ), incomingVelocity_,
                outgoingRotationAngle_, periapsisRadius_, swingbyDeltaV_ );

    totalNodeDeltaV_ = swingbyDeltaV_;

}


SwingbyWithFreeIncomingFixedOutgoingVelocity::SwingbyWithFreeIncomingFixedOutgoingVelocity(
           const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
           const double centralBodyGravitationalParameter,
           const std::function< Eigen::Vector3d( ) > outgoingVelocityFunction ):
           TransferNode( nodeEphemeris, swingby ),
       centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
       outgoingVelocityFunction_( outgoingVelocityFunction )
{ }

bool SwingbyWithFreeIncomingFixedOutgoingVelocity::nodeComputesOutgoingVelocity( )
{
   return false;
}

bool SwingbyWithFreeIncomingFixedOutgoingVelocity::nodeComputesIncomingVelocity( )
{
   return true;
}

void SwingbyWithFreeIncomingFixedOutgoingVelocity::computeNode( )
{
    if( nodeParameters_.rows( ) != 4 )
    {
        throw std::runtime_error( "Error when computing SwingbyWithFreeIncomingFixedOutgoingVelocity, incorrect input size" );
    }

    nodeTime_ = nodeParameters_( 0 );
    periapsisRadius_ = nodeParameters_( 1 );
    incomingRotationAngle_ = nodeParameters_( 2 );
    swingbyDeltaV_ = nodeParameters_( 3 );

    updateNodeState( nodeTime_ );

    outgoingVelocity_ = outgoingVelocityFunction_( );

    // Backward propagate the gravity assist
    incomingVelocity_ = mission_segments::calculatePoweredGravityAssistIncomingVelocity(
            centralBodyGravitationalParameter_, nodeState_.segment<3>(3), outgoingVelocity_,
            incomingRotationAngle_, periapsisRadius_, swingbyDeltaV_);

    totalNodeDeltaV_ = swingbyDeltaV_;
}

} // namespace mission_segments

} // namespace tudat


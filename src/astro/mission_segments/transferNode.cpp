#include "tudat/astro/mission_segments/transferNode.h"

namespace tudat
{
namespace mission_segments
{


TransferNode::TransferNode(
        const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
        const TransferNodeTypes nodeType,
        const Eigen::VectorXd nodeParameters ):
    nodeEphemeris_( nodeEphemeris ), nodeType_( nodeType ), nodeParameters_( nodeParameters ){ }

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





DepartureWithFixedOutgoingVelocityNode::DepartureWithFixedOutgoingVelocityNode(
        const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
        const Eigen::VectorXd nodeParameters,
        const double centralBodyGravitationalParameter,
        const double departureSemiMajorAxis,
        const double departureEccentricity,
        const std::function< Eigen::Vector3d( ) > outgoingVelocityFunction ):
    TransferNode( nodeEphemeris, escape_and_departure, nodeParameters ),
    centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
    departureSemiMajorAxis_( departureSemiMajorAxis ),
    departureEccentricity_( departureEccentricity ),
    outgoingVelocityFunction_( outgoingVelocityFunction )
{
//    updateNodeParameters( nodeParameters );
}

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

    outgoingVelocity_ = outgoingVelocityFunction_( );
    centralBodyVelocity_ = nodeEphemeris_->getCartesianState( nodeTime_ ).segment( 3, 3 );

    totalNodeDeltaV_ = mission_segments::computeEscapeOrCaptureDeltaV(
                centralBodyGravitationalParameter_,
                departureSemiMajorAxis_, departureEccentricity_,
                ( outgoingVelocity_ - centralBodyVelocity_ ).norm( ) );
}





DepartureWithFreeOutgoingVelocityNode::DepartureWithFreeOutgoingVelocityNode(
        const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
        const Eigen::VectorXd nodeParameters,
        const double centralBodyGravitationalParameter,
        const double departureSemiMajorAxis,
        const double departureEccentricity ):
    TransferNode( nodeEphemeris, escape_and_departure, nodeParameters ),
    centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
    departureSemiMajorAxis_( departureSemiMajorAxis ),
    departureEccentricity_( departureEccentricity )
{
//    updateNodeParameters( nodeParameters );
}

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

    centralBodyPosition_ = nodeEphemeris_->getCartesianState( nodeTime_ ).segment( 0, 3 );
    centralBodyVelocity_ = nodeEphemeris_->getCartesianState( nodeTime_ ).segment( 3, 3 );

    // Calculate unit vectors as described in [Vinko and Izzo, 2008].
    const Eigen::Vector3d unitVector1 = centralBodyVelocity_.normalized( );
    const Eigen::Vector3d unitVector3 = ( centralBodyPosition_.cross( centralBodyVelocity_ ) ).normalized( );
    const Eigen::Vector3d unitVector2 = unitVector3.cross( unitVector1 );


    // Calculate the velocity after departure as described in [Vinko and Izzo, 2008].
    outgoingVelocity_ = centralBodyVelocity_ +
            excessVelocityMagnitude_ * std::cos( excessVelocityInPlaneAngle_ ) *
            std::cos( excessVelocityOutOfPlaneAngle_ ) * unitVector1 +
            excessVelocityMagnitude_ * std::sin( excessVelocityInPlaneAngle_ ) *
            std::cos( excessVelocityOutOfPlaneAngle_ ) * unitVector2 +
            excessVelocityMagnitude_ * std::sin( excessVelocityOutOfPlaneAngle_ ) * unitVector3;

    totalNodeDeltaV_ = mission_segments::computeEscapeOrCaptureDeltaV(
                centralBodyGravitationalParameter_,
                departureSemiMajorAxis_, departureEccentricity_,
                ( outgoingVelocity_ - centralBodyVelocity_ ).norm( ) );
}


CaptureAndInsertionNode::CaptureAndInsertionNode(
        const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
        const Eigen::VectorXd nodeParameters,
        const double centralBodyGravitationalParameter,
        const double captureSemiMajorAxis,
        const double captureEccentricity,
        const std::function< Eigen::Vector3d( ) > incomingVelocityFunction ):
    TransferNode( nodeEphemeris, capture_and_insertion, nodeParameters ),
    centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
    captureSemiMajorAxis_( captureSemiMajorAxis ),
    captureEccentricity_( captureEccentricity ),
    incomingVelocityFunction_( incomingVelocityFunction )
{
//    updateNodeParameters( nodeParameters );
}

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

    incomingVelocity_ = incomingVelocityFunction_( );
    centralBodyVelocity_ = nodeEphemeris_->getCartesianState( nodeTime_ ).segment( 3, 3 );

    totalNodeDeltaV_ = mission_segments::computeEscapeOrCaptureDeltaV(
                centralBodyGravitationalParameter_,
                captureSemiMajorAxis_, captureEccentricity_,
                ( incomingVelocity_ - centralBodyVelocity_ ).norm( ) );
}


SwingbyWithFixedOutgoingVelocity::SwingbyWithFixedOutgoingVelocity(
        const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
        const Eigen::VectorXd nodeParameters,
        const double centralBodyGravitationalParameter,
        const double minimumPeriapsisRadius,
        const std::function< Eigen::Vector3d( ) > incomingVelocityFunction,
        const std::function< Eigen::Vector3d( ) > outgoingVelocityFunction ):
    TransferNode( nodeEphemeris, swingby, nodeParameters ),
    centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
    minimumPeriapsisRadius_( minimumPeriapsisRadius ),
    incomingVelocityFunction_( incomingVelocityFunction ),
    outgoingVelocityFunction_( outgoingVelocityFunction )
{
//    updateNodeParameters( nodeParameters );
}

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
    incomingVelocity_ = incomingVelocityFunction_( );
    outgoingVelocity_ = outgoingVelocityFunction_( );

    centralBodyVelocity_ = nodeEphemeris_->getCartesianState( nodeTime_ ).segment( 3, 3 );

    totalNodeDeltaV_ = calculateGravityAssistDeltaV(
                centralBodyGravitationalParameter_, centralBodyVelocity_,
                incomingVelocity_, outgoingVelocity_, minimumPeriapsisRadius_ );

}




SwingbyWithFreeOutgoingVelocity::SwingbyWithFreeOutgoingVelocity(
        const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
        const Eigen::VectorXd nodeParameters,
        const double centralBodyGravitationalParameter,
        const std::function< Eigen::Vector3d( ) > incomingVelocityFunction ):
    TransferNode( nodeEphemeris, swingby, nodeParameters ),
    centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
    incomingVelocityFunction_( incomingVelocityFunction )
{
//    updateNodeParameters( nodeParameters );
}

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

    incomingVelocity_ = incomingVelocityFunction_( );

    centralBodyVelocity_ = nodeEphemeris_->getCartesianState( nodeTime_ ).segment( 3, 3 );

    // Prepare the gravity assist propagator module.
    outgoingVelocity_ = mission_segments::calculatePoweredGravityAssistOutgoingVelocity(
                centralBodyGravitationalParameter_,
                centralBodyVelocity_, incomingVelocity_,
                outgoingRotationAngle_, periapsisRadius_, swingbyDeltaV_ );

    totalNodeDeltaV_ = swingbyDeltaV_;

}



} // namespace mission_segments

} // namespace tudat


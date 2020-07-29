#include "tudat/astro/mission_segments/createTransferTrajectory.h"

namespace tudat
{
namespace mission_segments
{

std::shared_ptr< TransferLegSettings > dsmVelocityBasedLeg( )
{
    return std::make_shared< TransferLegSettings >( dsm_velocity_based_leg );
}


std::shared_ptr< TransferLegSettings > dsmPositionBasedLeg( )
{
    return std::make_shared< TransferLegSettings >( dsm_position_based_leg );
}


std::shared_ptr< TransferLegSettings > unpoweredLeg( )
{
    return std::make_shared< TransferLegSettings >( unpowered_unperturbed_leg );
}

std::shared_ptr< TransferNodeSettings > escapeAndCaptureNode(
        const double departureSemiMajorAxis,
        const double departureEccentricity )
{
    return std::make_shared< EscapeAndDepartureNodeSettings >( departureSemiMajorAxis, departureEccentricity );
}

std::shared_ptr< TransferNodeSettings > swingbyNode(
        const double minimumPeriapsisDistance )
{
    return std::make_shared< SwingbyNodeSettings >( minimumPeriapsisDistance );
}

std::shared_ptr< CaptureAndInsertionNodeSettings > captureAndInsertionNode(
        const double captureSemiMajorAxis,
        const double captureEccentricity )
{
    return std::make_shared< CaptureAndInsertionNodeSettings >(
                captureSemiMajorAxis, captureEccentricity );
}

std::shared_ptr< TransferLeg > createTransferLeg(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< TransferLegSettings > legSettings,
        const std::string& departureBodyName,
        const std::string& arrivalBodyName,
        const std::string& centralBodyName,
        const double legStartTime,
        const double legEndTime,
        const Eigen::VectorXd& legFreeParameters,
        const std::shared_ptr< TransferNode > departureNode )
{
    if( bodyMap.count( centralBodyName ) == 0 )
    {
        throw std::runtime_error( "Error when making transfer leg, central body " + centralBodyName  + " not found." );
    }
    else if( bodyMap.at( centralBodyName )->getGravityFieldModel( ) == nullptr )
    {
        throw std::runtime_error( "Error when making transfer leg, central body " + centralBodyName  + " has no gravity field." );
    }

    if( bodyMap.count( departureBodyName ) == 0 )
    {
        throw std::runtime_error( "Error when making transfer leg, departure body " + departureBodyName  + " not found." );
    }
    else if( bodyMap.at( departureBodyName )->getEphemeris( ) == nullptr )
    {
        throw std::runtime_error( "Error when making transfer leg, departure body " + departureBodyName  + " has no ephemeris." );
    }

    if( bodyMap.count( arrivalBodyName ) == 0 )
    {
        throw std::runtime_error( "Error when making transfer leg, arrival body " + arrivalBodyName  + " not found." );
    }
    else if( bodyMap.at( arrivalBodyName )->getEphemeris( ) == nullptr )
    {
        throw std::runtime_error( "Error when making transfer leg, arrival body " + arrivalBodyName  + " has no ephemeris." );
    }


    const double centralBodyGravitationalParameter =
            bodyMap.at( centralBodyName )->getGravityFieldModel( )->getGravitationalParameter( );
    const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris =
            bodyMap.at( departureBodyName )->getEphemeris( );
    const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris =
            bodyMap.at( arrivalBodyName )->getEphemeris( );

    std::shared_ptr< TransferLeg >  transferLeg;
    switch( legSettings->legType_ )
    {
    case unpowered_unperturbed_leg:
    {
        transferLeg = std::make_shared< UnpoweredUnperturbedTransferLeg >(
                    departureBodyEphemeris, arrivalBodyEphemeris,
                    ( Eigen::VectorXd( 2 )<< legStartTime, legEndTime ).finished( ),
                    centralBodyGravitationalParameter );
        break;
    }
    case dsm_position_based_leg:
    {
        transferLeg = std::make_shared< DsmPositionBasedTransferLeg >(
                    departureBodyEphemeris, arrivalBodyEphemeris,
                    ( Eigen::VectorXd( 6 )<< legStartTime, legEndTime, legFreeParameters ).finished( ),
                    centralBodyGravitationalParameter );
        break;
    }
    case dsm_velocity_based_leg:
    {
        if( departureNode == nullptr )
        {
            throw std::runtime_error( "Error when making dsm_velocity_based_leg, no previous departure node" );
        }
        std::function< Eigen::Vector3d( ) > departureVelocityFunction =
                std::bind( &TransferNode::getOutgoingVelocity, departureNode );

        transferLeg = std::make_shared< DsmVelocityBasedTransferLeg >(
                    departureBodyEphemeris, arrivalBodyEphemeris,
                    ( Eigen::VectorXd( 3 )<< legStartTime, legEndTime, legFreeParameters ).finished( ),
                    centralBodyGravitationalParameter, departureVelocityFunction );
        break;
    }
    default:
        throw std::runtime_error( "Error when making transfer leg, leg type not recognized" );
    }
    return transferLeg;
}

std::shared_ptr< TransferNode > createTransferNode(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< TransferNodeSettings > nodeSettings,
        const std::string &nodeBodyName,
        const double nodeTime,
        const Eigen::VectorXd& nodeFreeParameters,
        const std::shared_ptr< TransferLeg > incomingTransferLeg,
        const std::shared_ptr< TransferLeg > outgoingTransferLeg,
        const bool nodeComputesOutgoingVelocity )
{

    if( bodyMap.count( nodeBodyName ) == 0 )
    {
        throw std::runtime_error( "Error when making transfer node, body " + nodeBodyName  + " not found." );
    }
    else if( bodyMap.at( nodeBodyName )->getGravityFieldModel( ) == nullptr )
    {
        throw std::runtime_error( "Error when making transfer node, body " + nodeBodyName  + " has no gravity field." );
    }
    else if( bodyMap.at( nodeBodyName )->getEphemeris( ) == nullptr )
    {
        throw std::runtime_error( "Error when making transfer node, body " + nodeBodyName  + " has no ephemeris." );
    }

    const double centralBodyGravitationalParameter =
            bodyMap.at( nodeBodyName )->getGravityFieldModel( )->getGravitationalParameter( );
    const std::shared_ptr< ephemerides::Ephemeris > centralBodyEphemeris =
            bodyMap.at( nodeBodyName )->getEphemeris( );

    std::shared_ptr< TransferNode > transferNode;
    switch( nodeSettings->nodeType_ )
    {
    case swingby:
    {
        if( nodeComputesOutgoingVelocity == 0 )
        {
            std::shared_ptr< SwingbyNodeSettings > swingbySettings =
                    std::dynamic_pointer_cast< SwingbyNodeSettings >( nodeSettings );
            if( swingbySettings == nullptr )
            {
                throw std::runtime_error( "Error when making swingby node, type is inconsistent" );
            }

            if( swingbySettings->minimumPeriapsisRadius_ != swingbySettings->minimumPeriapsisRadius_ )
            {
                throw std::runtime_error(
                            "Error when making swingby node, no minimum periapsis radius is provided" );
            }

            std::function< Eigen::Vector3d( ) > incomingVelocityFunction =
                    std::bind( &TransferLeg::getArrivalVelocity, incomingTransferLeg );
            std::function< Eigen::Vector3d( ) > outgoingVelocityFunction =
                    std::bind( &TransferLeg::getDepartureVelocity, outgoingTransferLeg );
            transferNode = std::make_shared< SwingbyWithFixedOutgoingVelocity >(
                        centralBodyEphemeris, ( Eigen::VectorXd( 1 )<< nodeTime ).finished( ),
                        centralBodyGravitationalParameter, swingbySettings->minimumPeriapsisRadius_,
                        incomingVelocityFunction, outgoingVelocityFunction );
        }
        else
        {
            std::function< Eigen::Vector3d( ) > incomingVelocityFunction =
                    std::bind( &TransferLeg::getArrivalVelocity, incomingTransferLeg );
            transferNode = std::make_shared< SwingbyWithFreeOutgoingVelocity >(
                        centralBodyEphemeris, ( Eigen::VectorXd( 4 )<< nodeTime, nodeFreeParameters ).finished( ),
                        centralBodyGravitationalParameter, incomingVelocityFunction );
        }


        break;
    }
    case escape_and_departure:
    {
        if( nodeComputesOutgoingVelocity == 0 )
        {
            std::shared_ptr< EscapeAndDepartureNodeSettings > escapeAndDepartureSettings =
                    std::dynamic_pointer_cast< EscapeAndDepartureNodeSettings >( nodeSettings );
            if( escapeAndDepartureSettings == nullptr )
            {
                throw std::runtime_error( "Error when making escape_and_departure node, type is inconsistent" );
            }


            std::function< Eigen::Vector3d( ) > outgoingVelocityFunction =
                    std::bind( &TransferLeg::getDepartureVelocity, outgoingTransferLeg );
            transferNode = std::make_shared< DepartureWithFixedOutgoingVelocityNode >(
                        centralBodyEphemeris, ( Eigen::VectorXd( 1 )<< nodeTime ).finished( ),
                        centralBodyGravitationalParameter,
                        escapeAndDepartureSettings->departureSemiMajorAxis_,
                        escapeAndDepartureSettings->departureEccentricity_,
                        outgoingVelocityFunction );
        }
        else
        {
            std::shared_ptr< EscapeAndDepartureNodeSettings > escapeAndDepartureSettings =
                    std::dynamic_pointer_cast< EscapeAndDepartureNodeSettings >( nodeSettings );
            if( escapeAndDepartureSettings == nullptr )
            {
                throw std::runtime_error( "Error when making escape_and_departure node, type is inconsistent" );
            }


            transferNode = std::make_shared< DepartureWithFreeOutgoingVelocityNode >(
                        centralBodyEphemeris, ( Eigen::VectorXd( 4 )<< nodeTime, nodeFreeParameters ).finished( ),
                        centralBodyGravitationalParameter,
                        escapeAndDepartureSettings->departureSemiMajorAxis_,
                        escapeAndDepartureSettings->departureEccentricity_ );
        }
    }
        break;
    case capture_and_insertion:
    {

        std::shared_ptr< CaptureAndInsertionNodeSettings > captureAndInsertionSettings =
                std::dynamic_pointer_cast< CaptureAndInsertionNodeSettings >( nodeSettings );
        if( captureAndInsertionSettings == nullptr )
        {
            throw std::runtime_error( "Error when making capture_and_insertion node, type is inconsistent" );
        }


        std::function< Eigen::Vector3d( ) > incomingVelocityFunction =
                std::bind( &TransferLeg::getArrivalVelocity, incomingTransferLeg );
        transferNode = std::make_shared< CaptureAndInsertionNode >(
                    centralBodyEphemeris, ( Eigen::VectorXd( 1 )<< nodeTime ).finished( ),
                    centralBodyGravitationalParameter,
                    captureAndInsertionSettings->captureSemiMajorAxis_,
                    captureAndInsertionSettings->captureEccentricity_,
                    incomingVelocityFunction );
        break;
    }
    default:
        throw std::runtime_error( "Error when making transfer node, node type not recognized" );
    }
    return transferNode;
}

std::shared_ptr< TransferTrajectory > createTransferTrajectory(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::vector< std::shared_ptr< TransferLegSettings > >& legSettings,
        const std::vector< std::shared_ptr< TransferNodeSettings > >& nodeSettings,
        const std::vector< std::string >& nodeIds,
        const std::string& centralBody,
        const std::vector< double >& nodeTimes,
        const std::vector< Eigen::VectorXd >& legFreeParameters,
        const std::vector< Eigen::VectorXd >& nodeFreeParameters )
{
    if( legSettings.size( ) + 1 != nodeSettings.size( ) )
    {
        throw std::runtime_error( "Error when making transfer trajectory, number of legs ( "
                                  + std::to_string( legSettings.size( ) ) +
                                  " ) and number of nodes ( "
                                  + std::to_string( nodeSettings.size( ) ) +
                                  " ) are incompatible" );
    }

    if( nodeIds.size( ) != nodeSettings.size( ) )
    {
        throw std::runtime_error( "Error when making transfer trajectory, number of nodes ( "
                                  + std::to_string( nodeSettings.size( ) ) +
                                  " ) and number of node names ( "
                                  + std::to_string( nodeIds.size( ) ) +
                                  " ) are incompatible" );
    }

    if( nodeIds.size( ) != nodeTimes.size( ) )
    {
        throw std::runtime_error( "Error when making transfer trajectory, number of node names ( "
                                  + std::to_string( nodeIds.size( ) ) +
                                  " ) and number of node times ( "
                                  + std::to_string( nodeTimes.size( ) ) +
                                  " ) are incompatible" );
    }




    std::vector< std::shared_ptr< TransferLeg > > legs;
    std::vector< std::shared_ptr< TransferNode > > nodes;


    for( unsigned int i = 0; i < legSettings.size( ); i++ )
    {
        if( legRequiresInputFromNode.at( legSettings.at( i )->legType_ ) )
        {
            nodes.push_back(
                        createTransferNode(
                            bodyMap, nodeSettings.at( i ), nodeIds.at( i ), nodeTimes.at( i ),
                            nodeFreeParameters.at( i ), ( i == 0 ? nullptr : legs.at( i -  1 ) ), nullptr, true ) );
            legs.push_back(
                        createTransferLeg(
                            bodyMap, legSettings.at( i ),
                            nodeIds.at( i ), nodeIds.at( i + 1 ), centralBody,
                            nodeTimes.at( i ), nodeTimes.at( i + 1 ),
                            legFreeParameters.at( i ), nodes.at( i ) ) );
        }
        else
        {
            legs.push_back(
                        createTransferLeg(
                            bodyMap, legSettings.at( i ),
                            nodeIds.at( i ), nodeIds.at( i + 1 ), centralBody,
                            nodeTimes.at( i ), nodeTimes.at( i + 1 ),
                            legFreeParameters.at( i ) ) );

            nodes.push_back(
                        createTransferNode(
                            bodyMap, nodeSettings.at( i ), nodeIds.at( i ), nodeTimes.at( i ),
                            nodeFreeParameters.at( i ), ( i == 0 ? nullptr : legs.at( i -  1 ) ), legs.at( i ), false ) );


        }
    }

    std::cout<<nodeSettings.size( )<<" "<<nodeFreeParameters.size( )<<std::endl;
    nodes.push_back(
                createTransferNode(
                    bodyMap, nodeSettings.at( legSettings.size( ) ),
                    nodeIds.at( legSettings.size( ) ),
                    nodeTimes.at( legSettings.size( ) ),
                    nodeFreeParameters.at( legSettings.size( ) ),
                    legs.at( legSettings.size( ) -  1 ), nullptr, false ) );

    return std::make_shared< TransferTrajectory >( legs, nodes );

}


void getParameterVectorDecompositionIndices(
        const std::vector< std::shared_ptr< TransferLegSettings > >& legSettings,
        const std::vector< std::shared_ptr< TransferNodeSettings > >& nodeSettings,
        std::vector< std::pair< int, int > >& legParameterIndices,
        std::vector< std::pair< int, int > >& nodeParameterIndices )
{
    // First N parameters are times/TOFs
    int currentParameterIndex = nodeSettings.size( );

    for( unsigned int i = 0; i < nodeSettings.size( ); i++ )
    {
        switch( nodeSettings.at( i )->nodeType_  )
        {
        case swingby:
            if( legRequiresInputFromNode.at( legSettings.at( i )->legType_ ) )
            {
                nodeParameterIndices.push_back( std::make_pair( currentParameterIndex, 3 ) );
                currentParameterIndex += 3;
            }
            else
            {
                nodeParameterIndices.push_back( std::make_pair( currentParameterIndex, 0 ) );
            }
            break;
        case escape_and_departure:
            if( legRequiresInputFromNode.at( legSettings.at( i )->legType_ ) )
            {
                nodeParameterIndices.push_back( std::make_pair( currentParameterIndex, 3 ) );
                currentParameterIndex += 3;
            }
            else
            {
                nodeParameterIndices.push_back( std::make_pair( currentParameterIndex, 0 ) );
            }
            break;
        case capture_and_insertion:
            nodeParameterIndices.push_back( std::make_pair( currentParameterIndex, 0 ) );
            break;
        }

        if( i != legSettings.size( ) )
        {
            switch( legSettings.at( i )->legType_  )
            {
            case unpowered_unperturbed_leg:
                legParameterIndices.push_back( std::make_pair( currentParameterIndex, 0 ) );
                break;
            case dsm_position_based_leg:
                legParameterIndices.push_back( std::make_pair( currentParameterIndex, 4 ) );
                currentParameterIndex += 4;
                break;
            case dsm_velocity_based_leg:
                legParameterIndices.push_back( std::make_pair( currentParameterIndex, 1 ) );
                currentParameterIndex += 1;
                break;
            }
        }
    }
}

void printTransferParameterDefinition(
        const std::vector< std::shared_ptr< TransferLegSettings > >& legSettings,
        const std::vector< std::shared_ptr< TransferNodeSettings > >& nodeSettings )
{
    std::vector< std::vector< std::string > > legParameterDefinitions;
    std::vector< std::vector< std::string > > nodeParameterDefinitions;

    for( unsigned int i = 0; i < legSettings.size( ); i++ )
    {
        std::vector< std::string  > currentLegIds;
        switch( legSettings.at( i )->legType_  )
        {
        case unpowered_unperturbed_leg:
            break;
        case dsm_position_based_leg:
            currentLegIds.push_back( "DSM (position-based) Time-of-flight fraction" );
            currentLegIds.push_back( "DSM (position-based) Dimensionless radius" );
            currentLegIds.push_back( "DSM (position-based) In-plane angle" );
            currentLegIds.push_back( "DSM (position-based) Out-of-plane angle" );
            break;
        case dsm_velocity_based_leg:
            currentLegIds.push_back( "DSM (velocity-based) Time-of-flight fraction" );
            break;
        }
        legParameterDefinitions.push_back( currentLegIds );
    }

    for( unsigned int i = 0; i < nodeSettings.size( ); i++ )
    {
        std::vector< std::string  > currentNodeIds;
        switch( nodeSettings.at( i )->nodeType_  )
        {
        case swingby:
            if( legRequiresInputFromNode.at( legSettings.at( i )->legType_ ) )
            {
                currentNodeIds.push_back( "Swingby periapsis" );
                currentNodeIds.push_back( "Swingby orbit-orientation rotation" );
                currentNodeIds.push_back( "Swingby Delta V" );
            }
            break;
        case escape_and_departure:
            if( legRequiresInputFromNode.at( legSettings.at( i )->legType_ ) )
            {
                currentNodeIds.push_back( "Excess velocity magnitude" );
                currentNodeIds.push_back( "Excess velocity in-plane angle" );
                currentNodeIds.push_back( "Excess velocity out-of-plane angle" );
            }
            break;
        case capture_and_insertion:
            break;
        }
        nodeParameterDefinitions.push_back( currentNodeIds );
    }

    int parameterIndex = 0;
    for( int i = 0; i < nodeSettings.size( ); i++ )
    {
        std::cout << "Parameter "<<parameterIndex<<": Node time "<<i<<std::endl;
        parameterIndex++;
    }

    for( unsigned int i = 0; i < legParameterDefinitions.size( ); i++ )
    {

        for( unsigned j = 0; j < nodeParameterDefinitions.at( i ).size( ); j++ )
        {
            std::cout << "Parameter "<<parameterIndex<<": Node "<<i<<" "<< nodeParameterDefinitions.at( i ).at( j )<<std::endl;
            parameterIndex++;
        }

        for( unsigned j = 0; j < legParameterDefinitions.at( i ).size( ); j++ )
        {
            std::cout << "Parameter "<<parameterIndex<<": Leg  "<<i<<" "<< legParameterDefinitions.at( i ).at( j )<<std::endl;
            parameterIndex++;
        }

    }
    std::cout<<std::endl;
}


} // namespace mission_segments

} // namespace tudat

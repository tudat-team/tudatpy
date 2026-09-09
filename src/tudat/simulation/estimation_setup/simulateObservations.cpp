#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"

namespace tudat
{

namespace simulation_setup
{

std::map< double, Eigen::VectorXd > getTargetAnglesAndRange( const simulation_setup::SystemOfBodies& bodies,
                                                             const std::pair< std::string, std::string > groundStationId,
                                                             const std::string& targetBody,
                                                             const std::vector< double > times,
                                                             const bool transmittingToTarget )
{
    using namespace observation_models;
    //    if( observingBody->getGroundStationMap( ).count( groundStationName ) == 0 )
    //    {
    //        throw std::runtime_error( "Error when computing elevating angle, station " + groundStationName +
    //                                  " not found on body " + observingBody->getBodyName( ) );
    //    }

    LinkEnds linkEnds;
    LinkEndType groundStationRole;
    if( transmittingToTarget )
    {
        linkEnds[ transmitter ] = groundStationId;
        linkEnds[ receiver ] = std::pair< std::string, std::string >( std::make_pair( targetBody, "" ) );
        groundStationRole = transmitter;
    }
    else
    {
        linkEnds[ receiver ] = groundStationId;
        linkEnds[ transmitter ] = std::pair< std::string, std::string >( std::make_pair( targetBody, "" ) );
        groundStationRole = receiver;
    }

    std::shared_ptr< ObservationModelSettings > oneWayRangeSettings =
            std::make_shared< ObservationModelSettings >( one_way_range, linkEnds );
    std::shared_ptr< ObservationSimulatorBase< double, double > > observationSimulator =
            createObservationSimulators( { oneWayRangeSettings }, bodies ).at( 0 );

    std::shared_ptr< TabulatedObservationSimulationSettings< double > > observationSimulationSettings =
            std::make_shared< TabulatedObservationSimulationSettings< double > >( one_way_range, linkEnds, times, groundStationRole );

    std::vector< std::shared_ptr< ObservationDependentVariableSettings > > dependentVariablesToSave;
    dependentVariablesToSave.push_back( std::make_shared< StationAngleObservationDependentVariableSettings >(
            station_elevation_angle, groundStationId, groundStationRole ) );
    dependentVariablesToSave.push_back( std::make_shared< StationAngleObservationDependentVariableSettings >(
            station_azimuth_angle, groundStationId, groundStationRole ) );
    addDependentVariablesToObservationSimulationSettings( { observationSimulationSettings }, dependentVariablesToSave, bodies );

    std::shared_ptr< observation_models::ObservationDataset< double, double > > observations =
            simulateObservationDataset( { observationSimulationSettings }, { observationSimulator }, bodies );

    const std::vector< double > observationTimes = observations->getObservationTimesForSet( 0 );
    const std::vector< Eigen::VectorXd > angles = observations->getDependentVariablesForSet( 0 );
    const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > observationValues = observations->getObservationsForSet( 0 );

    std::map< double, Eigen::VectorXd > anglesAndRange;
    for( unsigned int i = 0; i < observationTimes.size( ); ++i )
    {
        anglesAndRange[ observationTimes.at( i ) ] =
                ( Eigen::VectorXd( 3 ) << angles.at( i )( 0 ), angles.at( i )( 1 ), observationValues.at( i )( 0 ) ).finished( );
    }
    return anglesAndRange;
}

}  // namespace simulation_setup

}  // namespace tudat

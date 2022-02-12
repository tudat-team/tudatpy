#include "tudat/simulation/estimation_setup/simulateObservations.h"

namespace tudat
{

namespace simulation_setup
{


std::map< double, Eigen::VectorXd > getTargetAnglesAndRange(
        const simulation_setup::SystemOfBodies& bodies,
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
        linkEnds[ receiver ] = std::make_pair( targetBody, "" );
        groundStationRole = transmitter;
    }
    else
    {
        linkEnds[ receiver ] = groundStationId;
        linkEnds[ transmitter ] = std::make_pair( targetBody, "" );
        groundStationRole = receiver;
    }


    std::shared_ptr< ObservationModelSettings > oneWayRangeSettings =
            std::make_shared< ObservationModelSettings >( one_way_range, linkEnds );
    std::shared_ptr< ObservationSimulatorBase< double, double > > observationSimulator =
            createObservationSimulators( { oneWayRangeSettings }, bodies ).at( 0 );

    std::shared_ptr< TabulatedObservationSimulationSettings< double > > observationSimulationSettings =
            std::make_shared< TabulatedObservationSimulationSettings< double > >(
                one_way_range, linkEnds, times, groundStationRole );

    std::vector< std::shared_ptr< ObservationDependentVariableSettings > > dependentVariablesToSave;
    dependentVariablesToSave.push_back( std::make_shared< StationAngleObservationDependentVariableSettings >(
                                            station_elevation_angle, groundStationId, groundStationRole ) );
    dependentVariablesToSave.push_back( std::make_shared< StationAngleObservationDependentVariableSettings >(
                                            station_azimuth_angle, groundStationId, groundStationRole ) );
    addDependentVariablesToObservationSimulationSettings( { observationSimulationSettings }, dependentVariablesToSave, bodies );

    std::shared_ptr< observation_models::ObservationCollection< double, double > >  observations =
            simulateObservations( { observationSimulationSettings }, { observationSimulator}, bodies );

    std::shared_ptr< SingleObservationSet< double, double > > singleObservationSet =  observations->getSingleLinkAndTypeObservationSets(
            one_way_range, linkEnds ).at( 0 );
    std::map< double, Eigen::VectorXd > angles = singleObservationSet->getDependentVariableHistory( );
    std::map< double, Eigen::VectorXd > observationValues = singleObservationSet->getObservationsHistory( );

    std::map< double, Eigen::VectorXd > anglesAndRange;
    for( auto it: angles )
    {
        anglesAndRange[ it.first ] = ( Eigen::VectorXd( 3 )<<it.second( 0 ), it.second( 1 ), observationValues.at( it.first ) ).finished( );
    }
    return anglesAndRange;
}

}

}

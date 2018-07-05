/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include "Tudat/JsonInterface/UnitTests/unitTestSupport.h"
#include "Tudat/JsonInterface/Estimation/observation.h"
#include "Tudat/JsonInterface/Estimation/parameter.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"
#include "Tudat/SimulationSetup/EstimationSetup/createObservationModel.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_observation )

// Test 1: sphericalHarmonicGravity
BOOST_AUTO_TEST_CASE( test_json_acceleration_sphericalHarmonicGravity )
{
    using namespace tudat::simulation_setup;
    using namespace tudat::observation_models;
    using namespace tudat::estimatable_parameters;
    using namespace tudat::json_interface;

    // Specify initial and final time
    double initialEphemerisTime = 1.0E7;
    int numberOfSimulationDays = 30.0;
    double finalEphemerisTime = initialEphemerisTime + numberOfSimulationDays * 86400.0;
    double arcDuration = 3.01 * 86400.0;
    double arcOverlap = 3600.0;

    std::vector< double > arcStartTimes;
    double currentTime = initialEphemerisTime;
    while( currentTime <= finalEphemerisTime )
    {
        arcStartTimes.push_back( currentTime );
        currentTime += arcDuration;
    }

    // Create ground stations from geodetic positions.
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1" );
    groundStationNames.push_back( "Station2" );
    groundStationNames.push_back( "Station3" );


    // Create list of link ends in which station is receiver and in which station is transmitter (with spacecraft other link end).
    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;

    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ receiver ] = std::make_pair( "Vehicle", "" );
        stationTransmitterLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ transmitter ] = std::make_pair( "Vehicle", "" );
        stationReceiverLinkEnds.push_back( linkEnds );
    }

    // Define (arbitrarily) link ends to be used for 1-way range, 1-way doppler and angular position observables
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 1 ] );

    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 1 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 2 ] );

    linkEndsPerObservable[ angular_position ].push_back( stationReceiverLinkEnds[ 2 ] );
    linkEndsPerObservable[ angular_position ].push_back( stationTransmitterLinkEnds[ 1 ] );

    std::vector< boost::shared_ptr< EstimatableParameterSettings > > parameterNames;

    Eigen::VectorXd systemInitialState = Eigen::VectorXd::Zero( 6 * arcStartTimes.size( ) );
    parameterNames.push_back(
                boost::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
                    "Vehicle", systemInitialState, arcStartTimes, "Earth" ) );
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "global_metric", ppn_parameter_gamma ) );
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Vehicle", radiation_pressure_coefficient ) );
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Vehicle", constant_drag_coefficient ) );
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Earth", constant_rotation_rate ) );
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Earth", rotation_pole_position ) );
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Earth", ground_station_position, "Station1" ) );
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Earth", ground_station_position, "Station2" ) );
    parameterNames.push_back( boost::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                  linkEndsPerObservable.at( one_way_range ).at( 0 ), one_way_range, true ) );
    parameterNames.push_back( boost::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                  linkEndsPerObservable.at( one_way_range ).at( 1 ), one_way_range, true ) );
    parameterNames.push_back( boost::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  2, 0, 8, 8, "Earth", spherical_harmonics_cosine_coefficient_block ) );
    parameterNames.push_back( boost::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  2, 1, 8, 8, "Earth", spherical_harmonics_sine_coefficient_block ) );
    nlohmann::json jsonObject;
    std::vector< boost::shared_ptr< EstimatableParameterSettings > > parameterNames2;

    to_json( jsonObject, parameterNames );
    from_json( jsonObject, parameterNames2 );


    std::ofstream outputFile( "/home/dominic/Software/numericalAstrodynamicsTudatBundle/tudatBundle/tudat/Tudat/JsonInterface/UnitTests/outputParameter.json" );
    outputFile << jsonObject.dump( 1 );
    outputFile.close( );

    std::cout<<jsonObject<<std::endl;

//    std::string linkEndString = "transmitter:(Earth, Station1); receiver:(Vehicle, )";

//    LinkEnds test = boost::lexical_cast< LinkEnds >( linkEndString );
//    for( auto it = test.begin( ); it != test.end( );it++ )
//    {
//        std::cout<<it->first<<" "<<it->second.first<<" "<<it->second.second<<std::endl;
//    }

}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat

/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <string>

#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Mathematics/Statistics/basicStatistics.h"
#include "Tudat/Mathematics/Statistics/randomVariableGenerator.h"
#include "Tudat/SimulationSetup/tudatEstimationHeader.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::observation_models;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;
using namespace tudat::basic_astrodynamics;
using namespace tudat::coordinate_conversions;
using namespace tudat::statistics;

BOOST_AUTO_TEST_SUITE( test_observation_noise_models )

double ignoreInputeVariable( boost::function< double( ) > inputFreeFunction, const double dummyInput )
{
    return inputFreeFunction( );
}

BOOST_AUTO_TEST_CASE( testObservationNoiseModels )
{
    //Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Moon" );

    // Specify initial time
    double initialEphemerisTime = double( 1.0E7 );

    // Create bodies needed in simulation
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames );

    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Creatre ground stations: same position, but different representation
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1" );
    groundStationNames.push_back( "Station2" );
    groundStationNames.push_back( "Station3" );

    createGroundStation( bodyMap.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );
    createGroundStation( bodyMap.at( "Earth" ), "Station2", ( Eigen::Vector3d( ) << 0.0, -0.55, 2.0 ).finished( ), geodetic_position );
    createGroundStation( bodyMap.at( "Earth" ), "Station3", ( Eigen::Vector3d( ) << 0.0, 0.05, 4.0 ).finished( ), geodetic_position );

    // Define parameters.
    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;

    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ receiver ] = std::make_pair( "Moon", "" );
        stationTransmitterLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ transmitter ] = std::make_pair( "Moon", "" );
        stationReceiverLinkEnds.push_back( linkEnds );
    }

    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 1 ] );

    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 1 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 2 ] );

    linkEndsPerObservable[ angular_position ].push_back( stationReceiverLinkEnds[ 2 ] );
    linkEndsPerObservable[ angular_position ].push_back( stationTransmitterLinkEnds[ 1 ] );

    double rangeBias1 = 3.0;
    double rangeBias2 = -7.2;

    observation_models::ObservationSettingsMap observationSettingsMap;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;

        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            boost::shared_ptr< ObservationBiasSettings > biasSettings;
            if( ( currentObservable == one_way_range ) && ( i == 0 ) )
            {
                biasSettings = boost::make_shared< ConstantRelativeObservationBiasSettings >(
                            Eigen::Vector1d::Constant( rangeBias1 ) );
            }
            else if( ( currentObservable == one_way_range ) && ( i == 1 ) )
            {
                biasSettings = boost::make_shared< ConstantRelativeObservationBiasSettings >(
                            Eigen::Vector1d::Constant( rangeBias2 ) );
            }

            observationSettingsMap.insert(
                        std::make_pair( currentLinkEndsList.at( i ),
                                        boost::make_shared< ObservationSettings >(
                                            currentObservable, boost::shared_ptr< LightTimeCorrectionSettings >( ),
                                            biasSettings ) ) );
        }
    }

    std::map< ObservableType,
    boost::shared_ptr< ObservationSimulatorBase< double, double > > >  observationSimulators =
            createObservationSimulators( observationSettingsMap, bodyMap );

    std::vector< double > baseTimeList;
    double observationTimeStart = initialEphemerisTime + 1000.0;
    double observationInterval = 10.0;
    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 5000; j++ )
        {
            baseTimeList.push_back( observationTimeStart + static_cast< double >( i ) * 86400.0 +
                                    static_cast< double >( j ) * observationInterval );
        }
    }

    std::map< ObservableType, std::map< LinkEnds, boost::shared_ptr< ObservationSimulationTimeSettings< double > > > >
            measurementSimulationInput;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            measurementSimulationInput[ currentObservable ][ currentLinkEndsList.at( i ) ] =
                    boost::make_shared< TabulatedObservationSimulationTimeSettings< double > >(
                        receiver, baseTimeList );
        }
    }


    typedef Eigen::Matrix< double, Eigen::Dynamic, 1 > ObservationVectorType;
    typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< double >, LinkEndType > > >
            SingleObservablePodInputType;
    typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;

    // Simulate observations
    PodInputDataType idealObservationsAndTimes = simulateObservations< double, double >(
                measurementSimulationInput, observationSimulators );

    double constantOffset = 0.1;
    double constantStandardDeviation = 0.5;
    boost::function< double( ) > inputFreeNoiseFunction = createBoostContinuousRandomVariableGeneratorFunction(
                normal_boost_distribution, boost::assign::list_of( constantOffset )( constantStandardDeviation ), 0.0 );
    boost::function< double( const double ) > noiseFunction =
            boost::bind( &ignoreInputeVariable, inputFreeNoiseFunction, _1 );

    PodInputDataType constantNoiseObservationsAndTimes = simulateObservationsWithNoise< double, double >(
                measurementSimulationInput, observationSimulators, noiseFunction );

    std::map< ObservableType, std::map< LinkEnds, std::vector< double > > > observationDifference;
    std::map< ObservableType, std::map< LinkEnds, double > > meanObservationDifference;
    std::map< ObservableType, std::map< LinkEnds, double > > standardDeviationObservationDifference;

    for( PodInputDataType::const_iterator dataIterator = constantNoiseObservationsAndTimes.begin( );
         dataIterator != constantNoiseObservationsAndTimes.end( ); dataIterator++ )
    {
        for( SingleObservablePodInputType::const_iterator innerDataIterator = dataIterator->second.begin( );
            innerDataIterator != dataIterator->second.end( ); innerDataIterator++ )
        {
            Eigen::VectorXd dataDifference = innerDataIterator->second.first -
                    idealObservationsAndTimes.at( dataIterator->first ).at( innerDataIterator->first ).first;
            meanObservationDifference[ dataIterator->first ][ innerDataIterator->first ] =
                    computeAverageOfVectorComponents( dataDifference );
            standardDeviationObservationDifference[ dataIterator->first ][ innerDataIterator->first ] =
                    computeStandardDeviationOfVectorComponents( dataDifference );
            std::cout<<meanObservationDifference[ dataIterator->first ][ innerDataIterator->first ]<<" "<<
                                standardDeviationObservationDifference[ dataIterator->first ][ innerDataIterator->first ]<<std::endl;
        }
    }

    std::map< ObservableType, double > constantOffsets;
    constantOffsets[ one_way_range ] = -0.2;
    constantOffsets[ one_way_doppler ] = 2.4E-8;
    constantOffsets[ angular_position ] = 3.0E-7;

    std::map< ObservableType, double > constantStandardDeviations;
    constantOffsets[ one_way_range ] = 2.4;
    constantOffsets[ one_way_doppler ] = 7.5E-7;
    constantOffsets[ angular_position ] = -6.3E-6;

}

BOOST_AUTO_TEST_SUITE_END( )

}

}


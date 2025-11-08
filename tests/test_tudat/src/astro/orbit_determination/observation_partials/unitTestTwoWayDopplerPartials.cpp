/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <boost/lambda/lambda.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/astro/observation_models/oneWayDopplerObservationModel.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationRate.h"
#include "tudat/simulation/estimation_setup/createObservationPartials.h"
#include "tudat/support/numericalObservationPartial.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/math/basic/numericalDerivative.h"
#include "tudat/support/observationPartialTestFunctions.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_one_way_observation_partials )

//! Test partial derivatives of one-way doppler observable, using general test suite of observation partials.
BOOST_AUTO_TEST_CASE( testTwoWayDopplerPartials )
{
    using namespace tudat::gravitation;
    using namespace tudat::ground_stations;
    using namespace tudat::ephemerides;
    using namespace tudat::observation_models;
    using namespace tudat::simulation_setup;
    using namespace tudat::spice_interface;
    using namespace tudat::observation_partials;
    using namespace tudat::estimatable_parameters;

    // Define and create ground stations.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.resize( 2 );
    groundStations[ 0 ] = std::make_pair( "Earth", "Graz" );
    groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );

    // Test partials with constant ephemerides (allows test of position partials)
    for( unsigned int normalizeObservable = 0; normalizeObservable < 2; normalizeObservable++ )
    {
        for( unsigned int useFrequency = 0; useFrequency < 2; useFrequency++ )
        {
            std::cout<<"USE FREQUENCY "<<useFrequency<<std::endl;
            // Create environment
            SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true );

            std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings;

            if( useFrequency == 1 )
            {
                bodies.at( "Earth" )->getGroundStation( "Graz" )->getVehicleSystems( )->setTransponderTurnaroundRatio( );
                bodies.at( "Mars" )->getGroundStation( "MSL" )->setTransmittingFrequencyCalculator(
                        std::make_shared< ground_stations::ConstantFrequencyInterpolator >( 5.0E9 ) );
                ancilliarySettings = std::make_shared< observation_models::ObservationAncilliarySimulationSettings >( );
                ancilliarySettings->setAncilliaryDoubleVectorData(
                        frequency_bands, { static_cast< double >( x_band ), static_cast< double >( x_band ) } );
                ancilliarySettings->setAncilliaryDoubleData( doppler_reference_frequency, 0.0 );
                ancilliarySettings->setAncilliaryDoubleData( reception_reference_frequency_band,
                                                            convertFrequencyBandToDouble( x_band ) );
            }

            // Set link ends for observation model
            LinkDefinition linkEnds;
            linkEnds[ transmitter ] = groundStations[ 1 ];
            linkEnds[ reflector1 ] = groundStations[ 0 ];
            linkEnds[ receiver ] = groundStations[ 1 ];

            unsigned int caseStartIndex = ( useFrequency == 0 ) ? 0 : 1;
            for( unsigned int estimationCase = caseStartIndex; estimationCase < 2; estimationCase++ )
            {
                std::cout << "Case A " << normalizeObservable << " " << estimationCase << std::endl;
                // Generate one-way doppler model
                std::shared_ptr< ObservationModel< 1 > > twoWayDopplerModel;
                std::vector< std::string > perturbingBodies;
                perturbingBodies.push_back( "Earth" );
                if( estimationCase == 0 )
                {
                    twoWayDopplerModel = observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                            std::make_shared< observation_models::TwoWayDopplerObservationModelSettings >(
                                    linkEnds,
                                    std::vector< std::shared_ptr< LightTimeCorrectionSettings > >(
                                            { std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ) } ),
                                    nullptr,
                                    std::make_shared< LightTimeConvergenceCriteria >( ),
                                    static_cast< bool >( normalizeObservable ) ),
                            bodies );
                }
                else
                {
                    twoWayDopplerModel = observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                            std::make_shared< TwoWayDopplerObservationModelSettings >(
                                    std::make_shared< observation_models::OneWayDopplerObservationModelSettings >(
                                            getUplinkFromTwoWayLinkEnds( linkEnds ),
                                            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ),
                                            std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Mars" ),
                                            std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ),
                                            nullptr,
                                            std::make_shared< LightTimeConvergenceCriteria >( ),
                                            static_cast< bool >( normalizeObservable ) ),
                                    std::make_shared< observation_models::OneWayDopplerObservationModelSettings >(
                                            getDownlinkFromTwoWayLinkEnds( linkEnds ),
                                            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ),
                                            std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ),
                                            std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Mars" ),
                                            nullptr,
                                            std::make_shared< LightTimeConvergenceCriteria >( ),
                                            static_cast< bool >( normalizeObservable ) ),
                                    ( useFrequency == 0 ) ? two_way_doppler : doppler_measured_frequency ),
                            bodies );
                }

                // Create parameter objects.
                std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet;
                Eigen::VectorXd parameterPerturbationMultipliers = Eigen::Vector4d::Constant( 1.0 );
                fullEstimatableParameterSet = createEstimatableParameters( bodies, 1.1E7 );

                printEstimatableParameterEntries( fullEstimatableParameterSet );

                testObservationPartials< 1 >( twoWayDopplerModel,
                                              bodies,
                                              fullEstimatableParameterSet,
                                              linkEnds,
                                              ( useFrequency == 0 ) ? two_way_doppler : doppler_measured_frequency,
                                              ( useFrequency == 0 ) ? 1.0E-5 : 1.0E-4,
                                              true,
                                              true,
                                              10.0,
                                              parameterPerturbationMultipliers,
                                              ancilliarySettings,
                                              1.1E7,
                                              100.0 );
            }
        }

        // Test partials with real ephemerides (without test of position partials)
        for( unsigned int useFrequency = 0; useFrequency < 2; useFrequency++ )
        {
            std::cout<<"USE FREQUENCY "<<useFrequency<<std::endl;

            // Create environment
            SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, false );

            std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings;

            if( useFrequency == 1 )
            {
                bodies.at( "Earth" )->getGroundStation( "Graz" )->getVehicleSystems( )->setTransponderTurnaroundRatio( );
                bodies.at( "Mars" )->getGroundStation( "MSL" )->setTransmittingFrequencyCalculator(
                        std::make_shared< ground_stations::ConstantFrequencyInterpolator >( 5.0E9 ) );
                ancilliarySettings = std::make_shared< observation_models::ObservationAncilliarySimulationSettings >( );
                ancilliarySettings->setAncilliaryDoubleVectorData(
                        frequency_bands, { static_cast< double >( x_band ), static_cast< double >( x_band ) } );
                ancilliarySettings->setAncilliaryDoubleData( doppler_reference_frequency, 0.0 );
                ancilliarySettings->setAncilliaryDoubleData( reception_reference_frequency_band,
                                                             convertFrequencyBandToDouble( x_band ) );
            }

            // Set link ends for observation model
            LinkDefinition linkEnds;
            linkEnds[ transmitter ] = groundStations[ 1 ];
            linkEnds[ reflector1 ] = groundStations[ 0 ];
            linkEnds[ receiver ] = groundStations[ 1 ];

            unsigned int caseStartIndex = ( useFrequency == 0 ) ? 0 : 1;
            for( unsigned int estimationCase = caseStartIndex; estimationCase < 2; estimationCase++ )
            {
                std::cout << "Case B " << normalizeObservable << " " << estimationCase << std::endl;
                // Generate two-way doppler model
                std::shared_ptr< ObservationModel< 1 > > twoWayDopplerModel;
                std::vector< std::string > perturbingBodies;
                perturbingBodies.push_back( "Earth" );
                if( estimationCase == 0 )
                {
                    twoWayDopplerModel = observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                            std::make_shared< observation_models::TwoWayDopplerObservationModelSettings >(
                                    linkEnds,
                                    std::vector< std::shared_ptr< LightTimeCorrectionSettings > >(
                                            { std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ) } ),
                                    nullptr,
                                    std::make_shared< LightTimeConvergenceCriteria >( ),
                                    static_cast< bool >( normalizeObservable ) ),
                            bodies );
                }
                else
                {
                    twoWayDopplerModel = observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                            std::make_shared< TwoWayDopplerObservationModelSettings >(
                                    std::make_shared< OneWayDopplerObservationModelSettings >(
                                            getUplinkFromTwoWayLinkEnds( linkEnds ),
                                            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ),
                                            std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Mars" ),
                                            std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ) ),
                                    std::make_shared< OneWayDopplerObservationModelSettings >(
                                            getDownlinkFromTwoWayLinkEnds( linkEnds ),
                                            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ),
                                            std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ),
                                            std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Mars" ) ),
                                    ( useFrequency == 0 ) ? two_way_doppler : doppler_measured_frequency ),
                            bodies );
                }
                // Create parameter objects.
                std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet;
                Eigen::VectorXd parameterPerturbationMultipliers = Eigen::Vector4d::Constant( 1.0 );
                fullEstimatableParameterSet = createEstimatableParameters( bodies, 1.1E7 );

                printEstimatableParameterEntries( fullEstimatableParameterSet );

                testObservationPartials< 1 >( twoWayDopplerModel,
                                              bodies,
                                              fullEstimatableParameterSet,
                                              linkEnds,
                                              ( useFrequency == 0 ) ? two_way_doppler : doppler_measured_frequency,
                                              1.0E-4,
                                              false,
                                              true,
                                              1.0,
                                              parameterPerturbationMultipliers,
                                              ancilliarySettings,
                                              1.1E7,
                                              100.0 );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN



#include <limits>
#include <boost/test/unit_test.hpp>

#include <boost/lambda/lambda.hpp>
#include "tudat/basics/testMacros.h"

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/estimation_setup/createCartesianStatePartials.h"
#include "tudat/astro/orbit_determination/observation_partials/rotationMatrixPartial.h"
#include "tudat/astro/ephemerides/fullPlanetaryRotationModel.h"
#include "tudat/astro/orbit_determination/observation_partials/clockParameterPartials.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/estimation_setup/createClockPartials.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"
#include "tudat/io/basicInputOutput.h"


namespace tudat
{
namespace unit_tests
{

//Using declarations.
using namespace tudat::observation_models;
using namespace tudat::estimatable_parameters;
using namespace tudat::interpolators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::system_models;
using namespace tudat::observation_partials;
using namespace tudat::orbit_determination;
using namespace tudat::ground_stations;

BOOST_AUTO_TEST_SUITE( test_clock_partials )

Eigen::MatrixXd getGlobalExpectedClockPartials( const double timeIntoArc )
{
    Eigen::MatrixXd expectedGlobalSatelliteCoefficientOneWayPartials = Eigen::MatrixXd::Zero( 1, 3 );
    expectedGlobalSatelliteCoefficientOneWayPartials( 0 ) = std::pow( static_cast< double >( timeIntoArc ), 2 );
    expectedGlobalSatelliteCoefficientOneWayPartials( 1 ) = std::pow( static_cast< double >( timeIntoArc ), 0 );
    expectedGlobalSatelliteCoefficientOneWayPartials( 2 ) = std::pow( static_cast< double >( timeIntoArc ), 1 );

    return expectedGlobalSatelliteCoefficientOneWayPartials;
}

Eigen::MatrixXd getArcwiseExpectedClockPartials( const double timeIntoArc )
{
    Eigen::MatrixXd expectedArcWiseSatelliteCoefficientPartials = Eigen::MatrixXd::Zero( 1, 12 );
    expectedArcWiseSatelliteCoefficientPartials.block( 0, 6, 1, 3 ) = getGlobalExpectedClockPartials( timeIntoArc );

    return expectedArcWiseSatelliteCoefficientPartials;

}
BOOST_AUTO_TEST_CASE( test_ClockPartials )
{
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = 1.2E7;

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define celestial bodies in simulation.
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames );
    SystemOfBodies bodyMap = createSystemOfBodies( bodySettings );
    bodyMap.createEmptyBody( "LAGEOS" );
    std::shared_ptr< Body > lageos = bodyMap.at( "LAGEOS" );

    // Create dummy constant state of satellite.
    Eigen::Vector6d lageosKeplerianElements;
    lageosKeplerianElements[ semiMajorAxisIndex ] = 8000.0E3;
    lageosKeplerianElements[ eccentricityIndex ] = 0.0044;
    lageosKeplerianElements[ inclinationIndex ] = 109.89 * mathematical_constants::PI / 180.0;
    lageosKeplerianElements[ argumentOfPeriapsisIndex ] = 259.35 * mathematical_constants::PI / 180.0;
    lageosKeplerianElements[ longitudeOfAscendingNodeIndex ] = 31.56 * mathematical_constants::PI / 180.0;
    lageosKeplerianElements[ trueAnomalyIndex ] = 1.0;
    Eigen::Vector6d dummyLageosState = convertKeplerianToCartesianElements(
            lageosKeplerianElements, getBodyGravitationalParameter("Earth" ) );
    lageos->setEphemeris( std::make_shared< ConstantEphemeris >(
            dummyLageosState, "Earth", "ECLIPJ2000" ) );

    // Create ground stations
    std::map< std::pair< std::string, std::string >, Eigen::Vector3d > groundStationsToCreate;
    groundStationsToCreate[ std::make_pair( "Earth", "Graz" ) ] =
            ( Eigen::Vector3d( ) << 1.7E6, -6.2E6, 1.3E5 ).finished( );
    groundStationsToCreate[ std::make_pair( "Earth", "Station2" ) ] =
        ( Eigen::Vector3d( ) << 2.7E6, 3.2E6, 0.1E5 ).finished( );
    createGroundStations( bodyMap, groundStationsToCreate );

    std::shared_ptr< GroundStation > grazStation = bodyMap.at( "Earth" )->getGroundStation( "Graz" );

    // Create satellite time error arcs.
    std::vector< Time > satelliteClockErrorArcTimes;
    double satelliteSingleArcLength = 0.5E5;
    Time currentTime = Time( initialEphemerisTime );
    while( currentTime < Time( finalEphemerisTime ) )
    {
        satelliteClockErrorArcTimes.push_back( currentTime );
        currentTime += satelliteSingleArcLength;
    }
    satelliteClockErrorArcTimes.push_back( finalEphemerisTime + satelliteSingleArcLength );
    int numberOfArcs = satelliteClockErrorArcTimes.size( );

    // Set errors for all arcs of satellite.
    std::vector< double > lageosArcPolynomialErrors;
    lageosArcPolynomialErrors.push_back( 0.0 );
    lageosArcPolynomialErrors.push_back( 0.0 );
    lageosArcPolynomialErrors.push_back( 0.0 );

    // Create and set timing system of satellite.
    std::shared_ptr< TimingSystem > lageosTimingSystem = std::make_shared< TimingSystem >( satelliteClockErrorArcTimes, lageosArcPolynomialErrors );
    std::shared_ptr< VehicleSystems > lageosSystemHardware = std::make_shared< VehicleSystems >( );
    lageosSystemHardware->setTimingSystem( lageosTimingSystem );
    lageos->setVehicleSystems( lageosSystemHardware );

    // Create graz time error arcs.
    std::vector< Time > grazClockErrorArcTimes;
    double grazSingleArcLength = 1.2E5;
    currentTime = Time( initialEphemerisTime );
    while( currentTime < Time( finalEphemerisTime ) )
    {
        grazClockErrorArcTimes.push_back( currentTime );
        currentTime += grazSingleArcLength;
    }
    grazClockErrorArcTimes.push_back( finalEphemerisTime + grazSingleArcLength );

    // Set errors for all arcs of ground station.
    std::vector< double > grazArcPolynomialErrors;
    grazArcPolynomialErrors.push_back( 0.0 );
    grazArcPolynomialErrors.push_back( 0.0 );
    grazArcPolynomialErrors.push_back( 0.0 );

    // Create and set timing system of graz.
    std::shared_ptr< TimingSystem > grazTimingSystem = std::make_shared< TimingSystem >( grazClockErrorArcTimes, grazArcPolynomialErrors );
    grazStation->setTimingSystem( grazTimingSystem );

    // Set link ends for test model.
    LinkEnds oneWayUplinkLinkEnds;
    oneWayUplinkLinkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", "Graz" ) );
    oneWayUplinkLinkEnds[ receiver ] = LinkEndId( std::make_pair( "LAGEOS", "" ) );

    LinkEnds oneWayDownlinkLinkEnds;
    oneWayDownlinkLinkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", "Graz" ) );
    oneWayDownlinkLinkEnds[ transmitter ] = LinkEndId( std::make_pair( "LAGEOS", "" ) );

    LinkEnds twoWayLinkEnds;
    twoWayLinkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", "Graz" ) );
    twoWayLinkEnds[ retransmitter ] = LinkEndId( std::make_pair( "LAGEOS", "" ) );
    twoWayLinkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", "Graz" ) );

    LinkEnds threeWayLinkEnds;
    threeWayLinkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", "Graz" ) );
    threeWayLinkEnds[ retransmitter ] = LinkEndId( std::make_pair( "LAGEOS", "" ) );
    threeWayLinkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", "Station2" ) );

    LinkEnds threeWayLinkEndsInverse;
    threeWayLinkEndsInverse[ receiver ] = LinkEndId( std::make_pair( "Earth", "Graz" ) );
    threeWayLinkEndsInverse[ retransmitter ] = LinkEndId( std::make_pair( "LAGEOS", "" ) );
    threeWayLinkEndsInverse[ transmitter ] = LinkEndId( std::make_pair( "Earth", "Station2" ) );

    // Create range model.
    std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList;
    biasSettingsList.push_back( std::make_shared< TiminigSystemBiasSettings >( "LAGEOS", "" ) );
    biasSettingsList.push_back( std::make_shared< TiminigSystemBiasSettings >( "Earth", "Graz" ) );
    std::shared_ptr< ObservationModelSettings > uplinkRangeSettings = std::make_shared< ObservationModelSettings >(
        one_way_range, oneWayUplinkLinkEnds, nullptr, std::make_shared< MultipleObservationBiasSettings >( biasSettingsList ) );
    std::shared_ptr< ObservationModelSettings > downlinkRangeSettings = std::make_shared< ObservationModelSettings >(
        one_way_range, oneWayDownlinkLinkEnds, nullptr, std::make_shared< MultipleObservationBiasSettings >( biasSettingsList ) );
    std::shared_ptr< ObservationModelSettings > twoWayRangeSettings = std::make_shared< ObservationModelSettings >(
        n_way_range, twoWayLinkEnds, nullptr, std::make_shared< MultipleObservationBiasSettings >( biasSettingsList ) );
    std::shared_ptr< ObservationModelSettings > threeWayRangeSettings = std::make_shared< ObservationModelSettings >(
        n_way_range, threeWayLinkEnds, nullptr, std::make_shared< MultipleObservationBiasSettings >( biasSettingsList ) );
    std::shared_ptr< ObservationModelSettings > threeWayRangeSettingsInverse = std::make_shared< ObservationModelSettings >(
        n_way_range, threeWayLinkEndsInverse, nullptr, std::make_shared< MultipleObservationBiasSettings >( biasSettingsList ) );

    std::shared_ptr< ObservationModel< 1, double, double > > uplinkRangeModel = ObservationModelCreator< 1, double, double >::createObservationModel(
        uplinkRangeSettings, bodyMap );
    std::shared_ptr< ObservationModel< 1, double, double > > downlinkRangeModel = ObservationModelCreator< 1, double, double >::createObservationModel(
        downlinkRangeSettings, bodyMap );    
    std::shared_ptr< ObservationModel< 1, double, double > > twoWayRangeModel = ObservationModelCreator< 1, double, double >::createObservationModel(
        twoWayRangeSettings, bodyMap );
    std::shared_ptr< ObservationModel< 1, double, double > > threeWayRangeModel = ObservationModelCreator< 1, double, double >::createObservationModel(
        threeWayRangeSettings, bodyMap );
    std::shared_ptr< ObservationModel< 1, double, double > > threeWayInverseRangeModel = ObservationModelCreator< 1, double, double >::createObservationModel(
        threeWayRangeSettingsInverse, bodyMap );

    // Get range without timing errors.
    double testTime = 1.05E7 + 2.0E4;
//    Eigen::VectorXd uplinkRangeWithoutTimeError = uplinkRangeModel->computeObservations( testTime, receiver );
//
//    // Get times into arcs of timing errors of station and satellite.
//    std::pair< Time, int > timeIntoCurrentSatelliteArc = lageosTimingSystem->getTimeIntoCurrentArcAndArcIndex( testTime );
//    std::pair< Time, int > timeIntoCurrentStationArc = grazTimingSystem->getTimeIntoCurrentArcAndArcIndex(
//        testTime - uplinkRangeWithoutTimeError.x( ) / physical_constants::SPEED_OF_LIGHT );

    // Create parameters from global and arcwise satellite clock corrections.
    std::shared_ptr< GlobalPolynomialClockCorrections > globalSatellitePolynomialClockCorrections =
        std::make_shared< GlobalPolynomialClockCorrections >(
            lageosTimingSystem,  std::vector< int >( {2, 0, 1} ), "LAGEOS", "" );
    std::shared_ptr< MultiArcClockCorrections > arcWiseSatellitePolynomialClockCorrections =
        std::make_shared< MultiArcClockCorrections >(
            lageosTimingSystem,
            std::vector< int >( { 2, 0, 1 } ),
            std::vector< int >( { 0, 1, 10, 37 } ),
            "LAGEOS", "" );

    // Create parameters from global and arcwise station clock corrections.
    std::shared_ptr<GlobalPolynomialClockCorrections> globalStationPolynomialClockCorrections =
        std::make_shared<GlobalPolynomialClockCorrections>(
            grazTimingSystem, std::vector<int>( { 2, 0, 1 } ), "Earth", "Graz" );
    std::shared_ptr<MultiArcClockCorrections> arcWiseStationPolynomialClockCorrections =
        std::make_shared<MultiArcClockCorrections>(
            grazTimingSystem,
            std::vector<int>( { 2, 0, 1 } ),
            std::vector<int>( { 1, 0, 4, 16 } ),
            "Earth", "Graz" );

    // Create timing partials wrt satellite clock corrections, and test against expeted values
    {
        int currentSatelliteArc = 10;
        double timeIntoSatelliteArc = testTime - satelliteClockErrorArcTimes.at( currentSatelliteArc );
        Eigen::MatrixXd expectedGlobalSatelliteCoefficientOneWayPartials = getGlobalExpectedClockPartials( timeIntoSatelliteArc );
        Eigen::MatrixXd expectedArcWiseSatelliteCoefficientPartials = getArcwiseExpectedClockPartials( timeIntoSatelliteArc );

        for( unsigned int clockPartialTestCase = 0; clockPartialTestCase < 5; clockPartialTestCase++ )
        {
            std::shared_ptr< ObservationModel< 1, double, double > > observationModel = nullptr;
            std::vector< int > satelliteLinkEndIndices;
            if( clockPartialTestCase == 0 )
            {
                observationModel = uplinkRangeModel;
                satelliteLinkEndIndices = { 1 };
            }
            else if( clockPartialTestCase == 1 )
            {
                observationModel = downlinkRangeModel;
                satelliteLinkEndIndices = { 0 };
            }
            else if( clockPartialTestCase == 2 )
            {
                observationModel = twoWayRangeModel;
                satelliteLinkEndIndices = { 1, 2 };
            }
            else if( clockPartialTestCase == 3 )
            {
                observationModel = threeWayRangeModel;
                satelliteLinkEndIndices = { 1, 2 };
            }
            else if( clockPartialTestCase == 4 )
            {
                observationModel = threeWayInverseRangeModel;
                satelliteLinkEndIndices = { 1, 2 };
            }

            ObservableType observableType = observationModel->getObservableType( );
            LinkEnds testLinkEnds = observationModel->getLinkEnds( );
            int numberOfPartialObjects = satelliteLinkEndIndices.size( );

            {
                std::map< int, std::shared_ptr< TimingPartial > > timingPartialsWrtGlobalSatelliteClockParameters =
                        createTimingPartialWrtClockProperty< Eigen::VectorXd >(
                            testLinkEnds, observableType, globalSatellitePolynomialClockCorrections,
                            getClockInducedBiases( observationModel->getObservationBiasCalculator( ) ) );
                std::map< int, std::shared_ptr< TimingPartial > > timingPartialsWrtArcWiseSatelliteClockParameters =
                        createTimingPartialWrtClockProperty< Eigen::VectorXd >(
                            testLinkEnds, observableType, arcWiseSatellitePolynomialClockCorrections,
                            getClockInducedBiases( observationModel->getObservationBiasCalculator( ) ) );

                BOOST_CHECK_EQUAL( timingPartialsWrtGlobalSatelliteClockParameters.size( ), numberOfPartialObjects );
                BOOST_CHECK_EQUAL( timingPartialsWrtArcWiseSatelliteClockParameters.size( ), numberOfPartialObjects );

                int counter = 0;
                for( auto it : timingPartialsWrtGlobalSatelliteClockParameters )
                {
                    BOOST_CHECK_EQUAL( it.first, satelliteLinkEndIndices.at( counter ) );
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( it.second->getPartialOfClockErrorWrtParameter( testTime ),
                                                       expectedGlobalSatelliteCoefficientOneWayPartials, std::numeric_limits< double >::epsilon( ) );
                    counter++;
                }

                counter = 0;
                for( auto it : timingPartialsWrtArcWiseSatelliteClockParameters )
                {
                    BOOST_CHECK_EQUAL( it.first, satelliteLinkEndIndices.at( counter ) );
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( it.second->getPartialOfClockErrorWrtParameter( testTime ),
                                                       expectedArcWiseSatelliteCoefficientPartials, std::numeric_limits< double >::epsilon( ) )
                    counter++;
                }
            }

            {
                std::shared_ptr< ObservationPartial< 1 > > partialWrtGlobalSatelliteClockParameters =
                    ObservationPartialWrtClockCreator< Eigen::VectorXd, 1 >::createPartialWrtClockProperty(
                        testLinkEnds, observableType, globalSatellitePolynomialClockCorrections,
                        getClockInducedBiases( observationModel->getObservationBiasCalculator( ) ) );

                std::shared_ptr< ObservationPartial< 1 > > partialWrtArcwiseSatelliteClockParameters =
                    ObservationPartialWrtClockCreator< Eigen::VectorXd, 1 >::createPartialWrtClockProperty(
                        testLinkEnds, observableType, arcWiseSatellitePolynomialClockCorrections,
                        getClockInducedBiases( observationModel->getObservationBiasCalculator( ) ) );

                std::vector< double > linkEndTimes;
                std::vector< Eigen::Vector6d > linkEndStates;
                observationModel->computeObservationsWithLinkEndData( testTime, receiver, linkEndTimes, linkEndStates );

                std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > rangePartialWrtGlobalSatelliteClockParameters =
                    partialWrtGlobalSatelliteClockParameters->calculatePartial( linkEndStates, linkEndTimes );
                std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > rangePartialWrtArcwiseSatelliteClockParameters =
                    partialWrtArcwiseSatelliteClockParameters->calculatePartial( linkEndStates, linkEndTimes );


                BOOST_CHECK_EQUAL( rangePartialWrtGlobalSatelliteClockParameters.size( ), numberOfPartialObjects );
                BOOST_CHECK_EQUAL( rangePartialWrtArcwiseSatelliteClockParameters.size( ), numberOfPartialObjects );

                int counter = 0;
                for( unsigned int partialIndex = 0; partialIndex < rangePartialWrtGlobalSatelliteClockParameters.size( ); partialIndex++ )
                {
                    auto it = rangePartialWrtGlobalSatelliteClockParameters.at( partialIndex );
                    BOOST_CHECK_CLOSE_FRACTION( it.second, linkEndTimes.at( satelliteLinkEndIndices.at( counter ) ), std::numeric_limits< double >::epsilon( ) );
                    double timeIntoCurrentArc = linkEndTimes.at( satelliteLinkEndIndices.at( counter ) ) - satelliteClockErrorArcTimes.at( currentSatelliteArc );
                    Eigen::MatrixXd expectedPartials = getGlobalExpectedClockPartials( timeIntoCurrentArc ) *
                        ( ( satelliteLinkEndIndices.at( counter ) % 2 == 0 ) ? -1.0 : 1.0 );
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( it.first, ( expectedPartials * physical_constants::SPEED_OF_LIGHT ), std::numeric_limits< double >::epsilon( ) );
                    counter++;
                }

                counter = 0;
                for( unsigned int partialIndex = 0; partialIndex < rangePartialWrtArcwiseSatelliteClockParameters.size( ); partialIndex++ )
                {
                    auto it = rangePartialWrtArcwiseSatelliteClockParameters.at( partialIndex );
                    BOOST_CHECK_CLOSE_FRACTION( it.second, linkEndTimes.at( satelliteLinkEndIndices.at( counter ) ), std::numeric_limits< double >::epsilon( ) );
                    double timeIntoCurrentArc = linkEndTimes.at( satelliteLinkEndIndices.at( counter ) ) - satelliteClockErrorArcTimes.at( currentSatelliteArc );
                    Eigen::MatrixXd expectedPartials = getArcwiseExpectedClockPartials( timeIntoCurrentArc ) *
                                                       ( ( satelliteLinkEndIndices.at( counter ) % 2 == 0 ) ? -1.0 : 1.0 );
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( it.first, ( expectedPartials * physical_constants::SPEED_OF_LIGHT ), std::numeric_limits< double >::epsilon( ) );
                    counter++;
                }
            }
        }
    }

    // Create timing partials wrt station clock corrections, and test against expected values
    {
        int currentStationArc = 4;
        double timeIntoStationArc = testTime - grazClockErrorArcTimes.at( currentStationArc );
        Eigen::MatrixXd expectedGlobalStationCoefficientOneWayPartials = getGlobalExpectedClockPartials( timeIntoStationArc );
        Eigen::MatrixXd expectedArcWiseStationCoefficientPartials = getArcwiseExpectedClockPartials( timeIntoStationArc );
        for( unsigned int clockPartialTestCase = 0; clockPartialTestCase < 5; clockPartialTestCase++ )
        {
            std::shared_ptr< ObservationModel< 1, double, double > > observationModel = nullptr;
            std::vector< int > stationLinkEndIndices;
            if( clockPartialTestCase == 0 )
            {
                observationModel = uplinkRangeModel;
                stationLinkEndIndices = { 0 };
            }
            else if( clockPartialTestCase == 1 )
            {
                observationModel = downlinkRangeModel;
                stationLinkEndIndices = { 1 };
            }
            else if( clockPartialTestCase == 2 )
            {
                observationModel = twoWayRangeModel;
                stationLinkEndIndices = { 0, 3 };
            }
            else if( clockPartialTestCase == 3 )
            {
                observationModel = threeWayRangeModel;
                stationLinkEndIndices = { 0 };
            }
            else if( clockPartialTestCase == 4 )
            {
                observationModel = threeWayInverseRangeModel;
                stationLinkEndIndices = { 3 };
            }

            ObservableType observableType = observationModel->getObservableType( );
            LinkEnds testLinkEnds = observationModel->getLinkEnds( );
            int numberOfPartialObjects = stationLinkEndIndices.size( );

            {
                std::map< int, std::shared_ptr< TimingPartial > > timingPartialsWrtGlobalStationClockParameters =
                    createTimingPartialWrtClockProperty< Eigen::VectorXd >(
                        testLinkEnds, observableType, globalStationPolynomialClockCorrections,
                        getClockInducedBiases( observationModel->getObservationBiasCalculator( ) ) );
                std::map< int, std::shared_ptr< TimingPartial > > timingPartialsWrtArcWiseStationClockParameters =
                    createTimingPartialWrtClockProperty< Eigen::VectorXd >(
                        testLinkEnds, observableType, arcWiseStationPolynomialClockCorrections,
                        getClockInducedBiases( observationModel->getObservationBiasCalculator( ) ) );

                BOOST_CHECK_EQUAL( timingPartialsWrtGlobalStationClockParameters.size( ), numberOfPartialObjects );
                BOOST_CHECK_EQUAL( timingPartialsWrtArcWiseStationClockParameters.size( ), numberOfPartialObjects );

                int counter = 0;
                for( auto it : timingPartialsWrtGlobalStationClockParameters )
                {
                    BOOST_CHECK_EQUAL( it.first, stationLinkEndIndices.at( counter ) );
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( it.second->getPartialOfClockErrorWrtParameter( testTime ),
                                                       expectedGlobalStationCoefficientOneWayPartials, std::numeric_limits< double >::epsilon( ) );
                    counter++;
                }

                counter = 0;
                for( auto it : timingPartialsWrtArcWiseStationClockParameters )
                {
                    BOOST_CHECK_EQUAL( it.first, stationLinkEndIndices.at( counter ) );
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( it.second->getPartialOfClockErrorWrtParameter( testTime ),
                                                       expectedArcWiseStationCoefficientPartials, std::numeric_limits< double >::epsilon( ) )
                    counter++;
                }
            }

            {
                std::shared_ptr< ObservationPartial< 1 > > partialWrtGlobalStationClockParameters =
                    ObservationPartialWrtClockCreator< Eigen::VectorXd, 1 >::createPartialWrtClockProperty(
                        testLinkEnds, observableType, globalStationPolynomialClockCorrections,
                        getClockInducedBiases( observationModel->getObservationBiasCalculator( ) ) );

                std::shared_ptr< ObservationPartial< 1 > > partialWrtArcwiseStationClockParameters =
                    ObservationPartialWrtClockCreator< Eigen::VectorXd, 1 >::createPartialWrtClockProperty(
                        testLinkEnds, observableType, arcWiseStationPolynomialClockCorrections,
                        getClockInducedBiases( observationModel->getObservationBiasCalculator( ) ) );

                std::vector< double > linkEndTimes;
                std::vector< Eigen::Vector6d > linkEndStates;
                observationModel->computeObservationsWithLinkEndData( testTime, receiver, linkEndTimes, linkEndStates );

                std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > rangePartialWrtGlobalStationClockParameters =
                    partialWrtGlobalStationClockParameters->calculatePartial( linkEndStates, linkEndTimes );
                std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > rangePartialWrtArcwiseStationClockParameters =
                    partialWrtArcwiseStationClockParameters->calculatePartial( linkEndStates, linkEndTimes );


                BOOST_CHECK_EQUAL( rangePartialWrtGlobalStationClockParameters.size( ), numberOfPartialObjects );
                BOOST_CHECK_EQUAL( rangePartialWrtArcwiseStationClockParameters.size( ), numberOfPartialObjects );

                int counter = 0;
                for( unsigned int partialIndex = 0; partialIndex < rangePartialWrtGlobalStationClockParameters.size( ); partialIndex++ )
                {
                    auto it = rangePartialWrtGlobalStationClockParameters.at( partialIndex );
                    BOOST_CHECK_CLOSE_FRACTION( it.second, linkEndTimes.at( stationLinkEndIndices.at( counter ) ), std::numeric_limits< double >::epsilon( ) );
                    double timeIntoCurrentArc = linkEndTimes.at( stationLinkEndIndices.at( counter ) ) - grazClockErrorArcTimes.at( currentStationArc );
                    Eigen::MatrixXd expectedPartials = getGlobalExpectedClockPartials( timeIntoCurrentArc ) *
                                                       ( ( stationLinkEndIndices.at( counter ) % 2 == 0 ) ? -1.0 : 1.0 );
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( it.first, ( expectedPartials * physical_constants::SPEED_OF_LIGHT ), std::numeric_limits< double >::epsilon( ) );
                    counter++;
                }

                counter = 0;
                for( unsigned int partialIndex = 0; partialIndex < rangePartialWrtArcwiseStationClockParameters.size( ); partialIndex++ )
                {
                    auto it = rangePartialWrtArcwiseStationClockParameters.at( partialIndex );
                    BOOST_CHECK_CLOSE_FRACTION( it.second, linkEndTimes.at( stationLinkEndIndices.at( counter ) ), std::numeric_limits< double >::epsilon( ) );
                    double timeIntoCurrentArc = linkEndTimes.at( stationLinkEndIndices.at( counter ) ) - grazClockErrorArcTimes.at( currentStationArc );
                    Eigen::MatrixXd expectedPartials = getArcwiseExpectedClockPartials( timeIntoCurrentArc ) *
                                                       ( ( stationLinkEndIndices.at( counter ) % 2 == 0 ) ? -1.0 : 1.0 );
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( it.first, ( expectedPartials * physical_constants::SPEED_OF_LIGHT ), std::numeric_limits< double >::epsilon( ) );
                    counter++;
                }
            }
        }
    }

}
BOOST_AUTO_TEST_SUITE_END( )

}

}



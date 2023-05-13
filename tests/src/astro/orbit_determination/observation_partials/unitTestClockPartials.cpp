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
    LinkEnds linkEndList;
    linkEndList[ transmitter ] = LinkEndId( std::make_pair( "Earth", "Graz" ) );
    linkEndList[ receiver ] = LinkEndId( std::make_pair( "LAGEOS", "" ) );

    // Link timing systems to observation biases
    std::vector< LinkEnds > allLinkEndsList;
    allLinkEndsList.push_back( linkEndList );
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_range ] = allLinkEndsList;

    // Create range model.
    std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList;
    biasSettingsList.push_back( std::make_shared< TiminigSystemBiasSettings >( "LAGEOS", "" ) );
    biasSettingsList.push_back( std::make_shared< TiminigSystemBiasSettings >( "Earth", "Graz" ) );
    std::shared_ptr< ObservationModelSettings > rangeSettings = std::make_shared< ObservationModelSettings >(
            one_way_range, linkEndList, nullptr, std::make_shared< MultipleObservationBiasSettings >( biasSettingsList ) );

    std::shared_ptr< ObservationModel< 1, double, double > > rangeModel = ObservationModelCreator< 1, double, double >::createObservationModel(
            rangeSettings, bodyMap );

    // Get range without timing errors.
    double testTime = 1.05E7 + 2.0E4;
    Eigen::VectorXd rangeWithoutTimeError = rangeModel->computeObservations( testTime, receiver );

    // Get times into arcs of timing errors of station and satellite.
    std::pair< Time, int > timeIntoCurrentSatelliteArc = lageosTimingSystem->getTimeIntoCurrentArcAndArcIndex(
                testTime );
    std::pair< Time, int > timeIntoCurrentStationArc = grazTimingSystem->getTimeIntoCurrentArcAndArcIndex(
                testTime - rangeWithoutTimeError.x( ) / physical_constants::SPEED_OF_LIGHT );


    // Create parameters from global and arcwise satellite clock corrections.
    std::shared_ptr< GlobalPolynomialClockCorrections > globalSatellitePolynomialClockCorrections =
            std::make_shared< GlobalPolynomialClockCorrections >( lageosTimingSystem,
                                                                  std::vector< int >( {2, 0, 1} ),
                                                                    "LAGEOS", "" );
    std::shared_ptr< MultiArcClockCorrections > arcWiseSatellitePolynomialClockCorrections =
            std::make_shared< MultiArcClockCorrections >(
                lageosTimingSystem,
                std::vector< int >( {2, 0, 1} ),
                std::vector< int >( { 1, 0, timeIntoCurrentSatelliteArc.second, static_cast< int >( satelliteClockErrorArcTimes.size( ) ) - 2 } ),
                "LAGEOS", "" );

    Eigen::MatrixXd expectedGlobalSatelliteCoefficientPartials = Eigen::MatrixXd::Zero( 1, 3 );
    expectedGlobalSatelliteCoefficientPartials( 0 ) = std::pow( static_cast< double >( timeIntoCurrentSatelliteArc.first ), 2 );
    expectedGlobalSatelliteCoefficientPartials( 1 ) = std::pow( static_cast< double >( timeIntoCurrentSatelliteArc.first ), 0 );
    expectedGlobalSatelliteCoefficientPartials( 2 ) = std::pow( static_cast< double >( timeIntoCurrentSatelliteArc.first ), 1 );

    Eigen::MatrixXd expectedArcWiseSatelliteCoefficientPartials = Eigen::MatrixXd::Zero( 1, 12 );
    expectedArcWiseSatelliteCoefficientPartials.block( 0, 6, 1, 3 ) = expectedGlobalSatelliteCoefficientPartials;

    // Create timing partials wrt satellite clock corrections.
    std::map< int, std::shared_ptr< TimingPartial > > timingPartialsWrtGlobalSatelliteClockParameters =
            createTimingPartialWrtClockProperty< Eigen::VectorXd >(
                linkEndList, one_way_range, globalSatellitePolynomialClockCorrections,
                getClockInducedBiases( rangeModel->getObservationBiasCalculator( ) ) );
    std::map< int, std::shared_ptr< TimingPartial > > timingPartialsWrtArcWiseSatelliteClockParameters =
            createTimingPartialWrtClockProperty< Eigen::VectorXd >(
                linkEndList, one_way_range, arcWiseSatellitePolynomialClockCorrections,
                getClockInducedBiases( rangeModel->getObservationBiasCalculator( ) ));

    BOOST_CHECK_EQUAL( timingPartialsWrtGlobalSatelliteClockParameters.size( ), 1 );
    BOOST_CHECK_EQUAL( timingPartialsWrtArcWiseSatelliteClockParameters.size( ), 1 );

    BOOST_CHECK_EQUAL( timingPartialsWrtGlobalSatelliteClockParameters.begin( )->first, 1 );
    BOOST_CHECK_EQUAL( timingPartialsWrtArcWiseSatelliteClockParameters.begin( )->first, 1 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( timingPartialsWrtGlobalSatelliteClockParameters.begin( )->second->getPartialOfClockErrorWrtParameter( testTime ),
                                       expectedGlobalSatelliteCoefficientPartials, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( timingPartialsWrtArcWiseSatelliteClockParameters.begin( )->second->getPartialOfClockErrorWrtParameter( testTime ),
                                       expectedArcWiseSatelliteCoefficientPartials, std::numeric_limits< double >::epsilon( ) );

    // Create range partials wrt satellite clock corrections.
    std::shared_ptr< ObservationPartial< 1 > > oneWayRangePartialWrtGlobalSatelliteClockParameters =
            ObservationPartialWrtClockCreator< Eigen::VectorXd, 1 >::createPartialWrtClockProperty(
                linkEndList, one_way_range, globalSatellitePolynomialClockCorrections,
                getClockInducedBiases( rangeModel->getObservationBiasCalculator( ) ) );
    std::shared_ptr< ObservationPartial< 1 > > oneWayRangePartialWrtArcWiseSatelliteClockParameters =
            ObservationPartialWrtClockCreator< Eigen::VectorXd, 1 >::createPartialWrtClockProperty(
                linkEndList, one_way_range, arcWiseSatellitePolynomialClockCorrections,
                getClockInducedBiases( rangeModel->getObservationBiasCalculator( ) ) );

    std::vector< double > linkEndTimes;
    std::vector< Eigen::Vector6d > linkEndStates;
    rangeModel->computeObservationsWithLinkEndData( testTime, receiver, linkEndTimes, linkEndStates );
    rangeModel->computeObservationsWithLinkEndData( testTime, receiver, linkEndTimes, linkEndStates );

    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > rangePartialWrtGlobalSatelliteClockParameters =
            oneWayRangePartialWrtGlobalSatelliteClockParameters->calculatePartial( linkEndStates, linkEndTimes );
    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > rangePartialWrtArcWiseSatelliteClockParameters =
            oneWayRangePartialWrtArcWiseSatelliteClockParameters->calculatePartial( linkEndStates, linkEndTimes );

    BOOST_CHECK_EQUAL( rangePartialWrtGlobalSatelliteClockParameters.size( ), 1 );
    BOOST_CHECK_EQUAL( rangePartialWrtArcWiseSatelliteClockParameters.size( ), 1 );

    BOOST_CHECK_CLOSE_FRACTION( rangePartialWrtGlobalSatelliteClockParameters[ 0 ].second, testTime, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( rangePartialWrtArcWiseSatelliteClockParameters[ 0 ].second, testTime, std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( rangePartialWrtGlobalSatelliteClockParameters[ 0 ].first,
                                       ( 1.0 * expectedGlobalSatelliteCoefficientPartials * physical_constants::SPEED_OF_LIGHT ),
                                       std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( rangePartialWrtArcWiseSatelliteClockParameters[ 0 ].first,
                                       ( 1.0 * expectedArcWiseSatelliteCoefficientPartials * physical_constants::SPEED_OF_LIGHT ),
                                       std::numeric_limits< double >::epsilon( ) );

    // Create parameters from global and arcwise station clock corrections.
    std::shared_ptr< GlobalPolynomialClockCorrections > globalStationPolynomialClockCorrections =
            std::make_shared< GlobalPolynomialClockCorrections >( grazTimingSystem, std::vector< int >( {2, 0, 1} ),
                                                                    "Earth", "Graz" );
    std::shared_ptr< MultiArcClockCorrections > arcWiseStationPolynomialClockCorrections =
            std::make_shared< MultiArcClockCorrections >(
                grazTimingSystem, std::vector< int >( {2, 0, 1} ),
                std::vector< int >( { 1, 0, timeIntoCurrentStationArc.second, static_cast< int >( grazClockErrorArcTimes.size( ) ) - 2 } ),
                "Earth", "Graz" );


    Eigen::MatrixXd expectedGlobalStationCoefficientPartials = Eigen::MatrixXd::Zero( 1, 3 );
    expectedGlobalStationCoefficientPartials( 0 ) = std::pow( static_cast< double >( timeIntoCurrentStationArc.first ), 2.0 );
    expectedGlobalStationCoefficientPartials( 1 ) = std::pow( static_cast< double >( timeIntoCurrentStationArc.first ), 0.0 );
    expectedGlobalStationCoefficientPartials( 2 ) = std::pow( static_cast< double >( timeIntoCurrentStationArc.first ), 1.0 );


    Eigen::MatrixXd expectedArcWiseStationCoefficientPartials = Eigen::MatrixXd::Zero( 1, 12 );
    expectedArcWiseStationCoefficientPartials.block( 0, 6, 1, 3 ) = expectedGlobalStationCoefficientPartials;

    // Create timing partials wrt station clock corrections.
    std::map< int, std::shared_ptr< TimingPartial > > timingPartialsWrtGlobalStationClockParameters =
            createTimingPartialWrtClockProperty< Eigen::VectorXd >(
                linkEndList, one_way_range, globalStationPolynomialClockCorrections,
                getClockInducedBiases( rangeModel->getObservationBiasCalculator( ) ) );
    std::map< int, std::shared_ptr< TimingPartial > > timingPartialsWrtArcWiseStationClockParameters =
            createTimingPartialWrtClockProperty< Eigen::VectorXd >(
                linkEndList, one_way_range, arcWiseStationPolynomialClockCorrections,
                getClockInducedBiases( rangeModel->getObservationBiasCalculator( ) ) );

    BOOST_CHECK_EQUAL( timingPartialsWrtGlobalStationClockParameters.size( ), 1 );
    BOOST_CHECK_EQUAL( timingPartialsWrtArcWiseStationClockParameters.size( ), 1 );

    BOOST_CHECK_EQUAL( timingPartialsWrtGlobalStationClockParameters.begin( )->first, 0 );
    BOOST_CHECK_EQUAL( timingPartialsWrtArcWiseStationClockParameters.begin( )->first, 0 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( timingPartialsWrtGlobalStationClockParameters.begin( )->second->getPartialOfClockErrorWrtParameter(
                                           testTime - rangeWithoutTimeError.x( ) / physical_constants::SPEED_OF_LIGHT ),
                                       expectedGlobalStationCoefficientPartials, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( timingPartialsWrtArcWiseStationClockParameters.begin( )->second->getPartialOfClockErrorWrtParameter(
                                           testTime - rangeWithoutTimeError.x( ) / physical_constants::SPEED_OF_LIGHT ),
                                       expectedArcWiseStationCoefficientPartials, std::numeric_limits< double >::epsilon( ) );

    // Create range partials wrt station clock corrections.
    std::shared_ptr< ObservationPartial< 1 > > oneWayRangePartialWrtGlobalStationClockParameters =
            ObservationPartialWrtClockCreator< Eigen::VectorXd, 1 >::createPartialWrtClockProperty(
                linkEndList, one_way_range, globalStationPolynomialClockCorrections,
                getClockInducedBiases( rangeModel->getObservationBiasCalculator( ) ) );
    std::shared_ptr< ObservationPartial< 1 > > oneWayRangePartialWrtArcWiseStationClockParameters =
            ObservationPartialWrtClockCreator< Eigen::VectorXd, 1 >::createPartialWrtClockProperty(
                linkEndList, one_way_range, arcWiseStationPolynomialClockCorrections,
                getClockInducedBiases( rangeModel->getObservationBiasCalculator( ) ) );


    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > rangePartialWrtGlobalStationClockParameters =
            oneWayRangePartialWrtGlobalStationClockParameters->calculatePartial( linkEndStates, linkEndTimes );
    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > rangePartialWrtArcWiseStationClockParameters =
            oneWayRangePartialWrtArcWiseStationClockParameters->calculatePartial( linkEndStates, linkEndTimes );

    BOOST_CHECK_EQUAL( rangePartialWrtGlobalStationClockParameters.size( ), 1 );
    BOOST_CHECK_EQUAL( rangePartialWrtArcWiseStationClockParameters.size( ), 1 );

    BOOST_CHECK_CLOSE_FRACTION( rangePartialWrtGlobalStationClockParameters[ 0 ].second, testTime
                                - rangeWithoutTimeError.x( ) / physical_constants::SPEED_OF_LIGHT, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( rangePartialWrtArcWiseStationClockParameters[ 0 ].second, testTime
                                - rangeWithoutTimeError.x( ) / physical_constants::SPEED_OF_LIGHT, std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( rangePartialWrtGlobalStationClockParameters[ 0 ].first,
                                       ( -expectedGlobalStationCoefficientPartials * physical_constants::SPEED_OF_LIGHT ),
                                       std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( rangePartialWrtArcWiseStationClockParameters[ 0 ].first,
                                       ( -expectedArcWiseStationCoefficientPartials * physical_constants::SPEED_OF_LIGHT ),
                                       std::numeric_limits< double >::epsilon( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

}

}



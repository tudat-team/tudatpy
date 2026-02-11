/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <limits>
#include <map>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <boost/filesystem.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/paths.hpp"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/astro/basic_astro/celestialBodyConstants.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/math/basic/leastSquaresEstimation.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/createMetric.h"
#include "tudat/simulation/environment_setup/createRelativisticTimeConverter.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/simulation.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::basic_astrodynamics;
using namespace tudat::input_output;
using namespace tudat::numerical_integrators;
using namespace tudat::simulation_setup;
using namespace tudat::unit_conversions;
using namespace tudat::propagators;

BOOST_AUTO_TEST_SUITE( test_iss_proper_time_rate )

BOOST_AUTO_TEST_CASE( ISS_proper_time_rate_iau )
{
    spice_interface::loadStandardSpiceKernels( );

    const std::string issCsvPath = "/Users/michael.plumaris/aces_data_analysis/Data/Relativistic/iss_tabulated.csv";
    Eigen::MatrixXd issData = input_output::readMatrixFromFile( issCsvPath, ",", "#" );
    const double initialEpoch = issData( 0, 0 );
    const double finalEpoch = initialEpoch + 10000;
    //const double finalEpoch   = issData( issData.rows( ) - 1, 0 );

    const double outputTimeStep = 10.0;
    const double ephemerisBuffer = physical_constants::JULIAN_DAY;

    std::vector< std::string > bodiesToCreate{ "Sun", "Earth", "Moon" };
    const auto globalFrameOrigin      = "Earth";
    const auto globalFrameOrientation = "J2000";
    auto bodySettings = getDefaultBodySettings(
                bodiesToCreate, initialEpoch - ephemerisBuffer, finalEpoch + ephemerisBuffer, globalFrameOrigin, globalFrameOrientation );

    // WGS84 Earth shape and high-accuracy GCRS->ITRS rotation
    const double flattening       = 1.0 / 298.257223563;
    const double equatorialRadius = 6378137.0;
    bodySettings.at( "Earth" )->shapeModelSettings =
            std::make_shared< simulation_setup::OblateSphericalBodyShapeSettings >( equatorialRadius, flattening );
    bodySettings.at( "Earth" )->rotationModelSettings =
            std::make_shared< simulation_setup::GcrsToItrsRotationModelSettings >( basic_astrodynamics::iau_2006, globalFrameOrientation );

    // Promote Earth gravity to degree/order 300 if available
    auto earthGravitySettings = std::dynamic_pointer_cast< simulation_setup::SphericalHarmonicsGravityFieldSettings >(
                bodySettings.at( "Earth" )->gravityFieldSettings );
    if( earthGravitySettings != nullptr )
    {
        const int targetDegree = 300;
        const int targetOrder  = 300;
        const int currentDegree = static_cast< int >( earthGravitySettings->getCosineCoefficients( ).rows( ) ) - 1;
        const int currentOrder  = static_cast< int >( earthGravitySettings->getCosineCoefficients( ).cols( ) ) - 1;
        const int degreeToCopy  = std::min( currentDegree, targetDegree );
        const int orderToCopy   = std::min( currentOrder, targetOrder );
        Eigen::MatrixXd cosine  = Eigen::MatrixXd::Zero( targetDegree + 1, targetOrder + 1 );
        Eigen::MatrixXd sine    = Eigen::MatrixXd::Zero( targetDegree + 1, targetOrder + 1 );
        cosine.block( 0, 0, degreeToCopy + 1, orderToCopy + 1 ) =
                earthGravitySettings->getCosineCoefficients( ).block( 0, 0, degreeToCopy + 1, orderToCopy + 1 );
        sine.block( 0, 0, degreeToCopy + 1, orderToCopy + 1 ) =
                earthGravitySettings->getSineCoefficients( ).block( 0, 0, degreeToCopy + 1, orderToCopy + 1 );
        earthGravitySettings = std::make_shared< simulation_setup::SphericalHarmonicsGravityFieldSettings >(
                earthGravitySettings->getGravitationalParameter( ),
                earthGravitySettings->getReferenceRadius( ),
                cosine, sine,
                earthGravitySettings->getAssociatedReferenceFrame( ) );
        earthGravitySettings->resetAssociatedReferenceFrame( "ITRS" );
        bodySettings.at( "Earth" )->gravityFieldSettings = earthGravitySettings;
    }

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // ISS tabulated ephemeris (epoch [s TDB], x y z vx vy vz in meters and m/s)
    std::map< double, Eigen::Vector6d > issStateHistory;
    for( int i = 0; i < issData.rows( ); ++i )
    {
        if( issData.cols( ) >= 7 )
        {
            Eigen::Vector6d state;
            state << issData( i, 1 ), issData( i, 2 ), issData( i, 3 ),
                     issData( i, 4 ), issData( i, 5 ), issData( i, 6 );
            issStateHistory[ issData( i, 0 ) ] = state;
        }
    }
    auto issInterpolator =
            std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::Vector6d > >( issStateHistory, 6 );
    auto issEphemeris = std::make_shared< ephemerides::TabulatedCartesianEphemeris< double, double > >(
            issInterpolator, globalFrameOrigin, globalFrameOrientation );
    bodies.createEmptyBody( "ISS" );
    bodies.getBody( "ISS" )->setEphemeris( issEphemeris );
    setGlobalFrameBodyEphemerides( bodies.getMap( ), globalFrameOrigin, globalFrameOrientation );

    // Initialize translational and rotational states at initial epoch
    for( const auto& bodyName : bodiesToCreate )
    {
        if( bodies.doesBodyExist( bodyName ) )
        {
            bodies.getBody( bodyName )->setStateFromEphemeris( initialEpoch );
        }
    }
    bodies.getBody( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialEpoch );
    bodies.getBody( "ISS" )->setStateFromEphemeris( initialEpoch );
//     {
//         auto earthBody = bodies.getBody( "Earth" );
//         auto earthGravity = std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
//                 earthBody->getGravityFieldModel( ) );
//         earthGravity->setRotationWrapper(
//                 std::make_shared< reference_frames::QuaternionRotationWrapper >(
//                         [earthBody]( )
//                         {
//                         try
//                         {
//                             const Eigen::Quaterniond q = Eigen::Quaterniond( earthBody->getCurrentRotationToLocalFrame( ) );
//                             std::cout << "[RotationWrapper] quat: [" << q.w() << ", " << q.x() << ", " << q.y() << ", " << q.z() << "]" << std::endl;
//                             return q;
//                         }
//                         catch( const std::exception& )
//                         {
//                             // Fallback to ephemeris at epoch zero if rotation not yet set
//                             earthBody->setCurrentRotationalStateToLocalFrameFromEphemeris( 0.0 );
//                             const Eigen::Quaterniond q = Eigen::Quaterniond( earthBody->getCurrentRotationToLocalFrame( ) );
//                             std::cout << "[RotationWrapper] quat (init): [" << q.w() << ", " << q.x() << ", " << q.y() << ", " << q.z() << "]" << std::endl;
//                             return q;
//                         }
//                         } ) );
//     }

    // Ground stations (height [m], latitude/longitude [deg]) - others commented for now
    std::map< std::string, Eigen::Vector3d > stationGeodetic;
    stationGeodetic[ "GT101" ] = ( Eigen::Vector3d( ) << unit_conversions::convertDegreesToRadians( 48.836 ),
                                 unit_conversions::convertDegreesToRadians( 2.3344 ), 137.5458 ).finished( );
    for( const auto& station : stationGeodetic )
    {
        createGroundStation(
                bodies.getBody( "Earth" ), station.first, station.second, coordinate_conversions::geodetic_position );
    }

    const double integratorStep = 10.0;
    auto integratorSettings = numerical_integrators::rungeKutta4Settings( initialEpoch, integratorStep );
    auto terminationSettings = std::make_shared< propagators::PropagationTimeTerminationSettings >( finalEpoch );

    // Chained: Earth second-order + GT101 topocentric, ISS second-order bodycentric
    const std::vector< std::string > topocentricPerturbingBodies{ "Sun", "Moon" };
    const Eigen::Matrix< double, Eigen::Dynamic, 1 > initialRelativisticState =
            Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( 1 );
    std::vector< std::shared_ptr< propagators::RelativisticTimeStatePropagatorSettings< double, double > > >
            bodyCentricToTopocentricConversionSettings;
    for( const auto& station : stationGeodetic )
    {
        bodyCentricToTopocentricConversionSettings.push_back(
                std::make_shared< propagators::BodycenteredToTopocentricTimePropagatorSettings< double, double > >(
                        std::make_pair( "Earth", station.first ),
                        false,
                        4,
                        true,
                        topocentricPerturbingBodies,
                        initialRelativisticState,
                        initialEpoch,
                        integratorSettings,
                        terminationSettings ) );
    }
    const std::vector< std::string > earthPerturbingBodies{ "Moon", "Sun" };
    const std::vector< std::string > issPerturbingBodies{ "Earth", "Sun", "Moon" };
    std::map< std::string, std::pair< int, int > > issSphericalHarmonics;
    issSphericalHarmonics[ "Earth" ] = std::make_pair( 300, 300 );
    std::map< std::string, std::shared_ptr< DirectRelativisticTimeConverterSettings<> > > converterSettings;
    converterSettings[ "Earth" ] = std::make_shared< DirectRelativisticTimeConverterSettings<> >(
            std::make_shared< propagators::SecondOrderBodyCenteredRelativisticTimeConverterSettings< double, double > >(
                    "Earth", earthPerturbingBodies, initialEpoch, integratorSettings, terminationSettings ),
            integratorSettings,
            bodyCentricToTopocentricConversionSettings );
    converterSettings[ "ISS" ] = std::make_shared< DirectRelativisticTimeConverterSettings<> >(
            std::make_shared< propagators::SecondOrderBodyCenteredRelativisticTimeConverterSettings< double, double > >(
                    "ISS", issPerturbingBodies, initialEpoch, integratorSettings, terminationSettings, issSphericalHarmonics ),
            integratorSettings );
    setRelativisticTimeConverters( bodies, converterSettings );

    auto earthTimeScaleConverter = bodies.getBody( "Earth" )->getTimeScaleConverter( );
    auto issTimeScaleConverter   = bodies.getBody( "ISS" )->getTimeScaleConverter( );
    BOOST_REQUIRE( earthTimeScaleConverter != nullptr );
    BOOST_REQUIRE( issTimeScaleConverter   != nullptr );

    for( const auto& station : stationGeodetic )
    {
        const std::string stationName = station.first;
        std::map< double, Eigen::VectorXd > conversionResults;
        for( double epoch = initialEpoch; epoch <= finalEpoch + std::numeric_limits< double >::epsilon( ); epoch += outputTimeStep )
        {
            Eigen::VectorXd current( 5 );
            const double tcbMinusStation = earthTimeScaleConverter->getTimeDifference(
                    barycentric_coordinate_time_scale, local_proper_time_scale, epoch, stationName );
            const double tcbMinusTcg = earthTimeScaleConverter->getTimeDifference(
                    barycentric_coordinate_time_scale, body_centered_coordinate_time_scale, epoch, stationName );
            const double tcgMinusStation = earthTimeScaleConverter->getTimeDifference(
                    body_centered_coordinate_time_scale, local_proper_time_scale, epoch, stationName );
            const double tcbMinusIss = issTimeScaleConverter->getTimeDifference(
                    barycentric_coordinate_time_scale, body_centered_coordinate_time_scale, epoch );
            const double ttMinusIss = ( tcbMinusIss - tcbMinusTcg ) +
                    physical_constants::LG_TIME_RATE_TERM * ( epoch - initialEpoch );

            current( 0 ) = tcbMinusTcg;                          // TCB - TCG
            current( 1 ) = tcgMinusStation;                      // TCG - Station
            current( 2 ) = tcbMinusIss;                          // TCB - ISS (body-centered)
            current( 3 ) = tcbMinusIss - tcbMinusStation;        // Station - ISS
            current( 4 ) = ttMinusIss;                           // TT - ISS
            conversionResults[ epoch ] = current;
        }

        const std::string outputDirectory = "/Users/michael.plumaris/aces_data_analysis/Data/Relativistic/";
        boost::filesystem::create_directories( outputDirectory );
        input_output::writeDataMapToTextFile(
                conversionResults,
                "test_ISS_proper_time_rate_iau_" + stationName + ".dat",
                outputDirectory, "", 16 );
    }

    BOOST_CHECK_EQUAL( 1, 1 );
}

BOOST_AUTO_TEST_CASE( ISS_proper_time_rate_metric )
{
    spice_interface::loadStandardSpiceKernels( );

    const std::string issCsvPath = "/Users/michael.plumaris/aces_data_analysis/Data/Relativistic/iss_tabulated.csv";
    Eigen::MatrixXd issData = input_output::readMatrixFromFile( issCsvPath, ",", "#" );
    const double initialEpoch = issData( 0, 0 );
    const double finalEpoch   = initialEpoch + 5000.0;
    const double ephemerisBuffer = physical_constants::JULIAN_DAY;

    std::vector< std::string > bodiesToCreate{ "Sun", "Earth", "Moon" };
    const auto globalFrameOrigin      = "Earth";
    const auto globalFrameOrientation = "J2000";
    auto bodySettings = getDefaultBodySettings(
                bodiesToCreate, initialEpoch - ephemerisBuffer, finalEpoch + ephemerisBuffer, globalFrameOrigin, globalFrameOrientation );

    bodySettings.at( "Earth" )->rotationModelSettings =
            std::make_shared< simulation_setup::GcrsToItrsRotationModelSettings >( basic_astrodynamics::iau_2006, globalFrameOrientation );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // ISS ephemeris from tabulated file (TDB)
    std::map< double, Eigen::Vector6d > issStateHistory;
    for( int i = 0; i < issData.rows( ); ++i )
    {
        if( issData.cols( ) >= 7 )
        {
            Eigen::Vector6d state;
            state << issData( i, 1 ), issData( i, 2 ), issData( i, 3 ),
                     issData( i, 4 ), issData( i, 5 ), issData( i, 6 );
            issStateHistory[ issData( i, 0 ) ] = state;
        }
    }
    auto issInterpolator =
            std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::Vector6d > >( issStateHistory, 6 );
    auto issEphemeris = std::make_shared< ephemerides::TabulatedCartesianEphemeris< double, double > >(
            issInterpolator, globalFrameOrigin, globalFrameOrientation );
    bodies.createEmptyBody( "ISS" );
    bodies.getBody( "ISS" )->setEphemeris( issEphemeris );
    setGlobalFrameBodyEphemerides( bodies.getMap( ), globalFrameOrigin, globalFrameOrientation );

    // Initialize states/rotation
    for( const auto& bodyName : bodiesToCreate )
    {
        if( bodies.doesBodyExist( bodyName ) )
        {
            bodies.getBody( bodyName )->setStateFromEphemeris( initialEpoch );
        }
    }
    bodies.getBody( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialEpoch );
    bodies.getBody( "ISS" )->setStateFromEphemeris( initialEpoch );

    // Metric settings: second-order Earth, first-order Sun/Moon.
    auto metricSettings = std::make_shared< SolarSystemSpaceTimeMetricSettings >(
            std::vector< std::string >{ "Sun", "Moon" },
            std::vector< std::string >{ "Earth" },
            std::map< std::string, std::pair< int, int > >( ),
            std::vector< std::string >( ),
            std::make_shared< relativity::PPNParameterSet >( 1.0, 1.0 ) );

    baseMetric = createSpaceTimeMetric( metricSettings, bodies );
    evaluatedMetricObjects.clear( );
    evaluatedMetricObjects[ std::make_pair( "ISS", "" ) ] = baseMetric->Clone( );

    const double integratorStep = 10.0;
    auto integratorSettings = numerical_integrators::rungeKutta4Settings( integratorStep );
    auto terminationSettings = std::make_shared< propagators::PropagationTimeTerminationSettings >( finalEpoch );

    auto directFromMetricSettings = std::make_shared< propagators::DirectRelativisticTimePropagatorSettings< double, double > >(
            std::make_pair( "ISS", "" ),
            initialEpoch,
            integratorSettings,
            terminationSettings );

    // Propagate ISS proper time directly from metric
    SingleArcDynamicsSimulator< double > dynamicsSimulator(
                bodies, directFromMetricSettings, true );
    const auto solution = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    BOOST_REQUIRE( !solution.empty( ) );

    const double coordinateDuration = finalEpoch - initialEpoch;
    const double properOffset       = solution.rbegin( )->second( 0 ); // integrates Δ(τ) relative to coordinate time
    const double rateOffset         = properOffset / coordinateDuration;

    // Proper time offset should be small (order 1e-8 of span) and finite
    BOOST_CHECK_SMALL( std::fabs( rateOffset ), 1.0E-6 );
    BOOST_CHECK( std::isfinite( properOffset ) );
}

BOOST_AUTO_TEST_CASE( geoid_clock_LG_rate )
{
    spice_interface::loadStandardSpiceKernels( );

    const double initialEpoch = 0.0;
    const double finalEpoch   = initialEpoch + 10.0 * physical_constants::JULIAN_DAY;
    const double ephemerisBuffer = physical_constants::JULIAN_DAY;

    const std::vector< std::string > bodiesToCreate{ "Sun", "Earth", "Moon" };
    const auto globalFrameOrigin      = "Earth";
    const auto globalFrameOrientation = "J2000";
    auto bodySettings = getDefaultBodySettings(
                bodiesToCreate, initialEpoch - ephemerisBuffer, finalEpoch + ephemerisBuffer, globalFrameOrigin, globalFrameOrientation );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodies.getMap( ), globalFrameOrigin, globalFrameOrientation );

    for( const auto& bodyName : bodiesToCreate )
    {
        if( bodies.doesBodyExist( bodyName ) )
        {
            bodies.getBody( bodyName )->setStateFromEphemeris( initialEpoch );
        }
    }
    bodies.getBody( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialEpoch );

    const std::string stationName = "GeoidStation";
    const Eigen::Vector3d stationGeodetic =
            ( Eigen::Vector3d( ) << 0.0, 0.0, 0.0 ).finished( ); // lat [rad], lon [rad], height [m]
    createGroundStation(
            bodies.getBody( "Earth" ), stationName, stationGeodetic, coordinate_conversions::geodetic_position );

    const double integratorStep = 60.0;
    auto integratorSettings = numerical_integrators::rungeKutta4Settings( integratorStep );
    auto terminationSettings = std::make_shared< propagators::PropagationTimeTerminationSettings >( finalEpoch );

    const std::vector< std::string > topocentricPerturbingBodies{ "Sun", "Moon" };
    const Eigen::Matrix< double, Eigen::Dynamic, 1 > initialRelativisticState =
            Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( 1 );
    std::vector< std::shared_ptr< propagators::RelativisticTimeStatePropagatorSettings< double, double > > >
            bodyCentricToTopocentricConversionSettings;
    bodyCentricToTopocentricConversionSettings.push_back(
            std::make_shared< propagators::BodycenteredToTopocentricTimePropagatorSettings< double, double > >(
                    std::make_pair( "Earth", stationName ),
                    false,
                    4,
                    true,
                    topocentricPerturbingBodies,
                    initialRelativisticState,
                    initialEpoch,
                    integratorSettings,
                    terminationSettings ) );

    const std::vector< std::string > earthPerturbingBodies{ "Moon", "Sun" };
    std::map< std::string, std::shared_ptr< DirectRelativisticTimeConverterSettings<> > > converterSettings;
    converterSettings[ "Earth" ] = std::make_shared< DirectRelativisticTimeConverterSettings<> >(
            std::make_shared< propagators::SecondOrderBodyCenteredRelativisticTimeConverterSettings< double, double > >(
                    "Earth", earthPerturbingBodies, initialEpoch, integratorSettings, terminationSettings ),
            integratorSettings,
            bodyCentricToTopocentricConversionSettings );
    setRelativisticTimeConverters( bodies, converterSettings );

    auto earthTimeScaleConverter = bodies.getBody( "Earth" )->getTimeScaleConverter( );
    BOOST_REQUIRE( earthTimeScaleConverter != nullptr );

    // Sample over 10 days with 60 s spacing and fit a slope of TCG-proper
    std::vector< double > timesVector;
    std::vector< double > resultDifferenceVector;
    for( double epoch = initialEpoch; epoch <= finalEpoch + std::numeric_limits< double >::epsilon( ); epoch += integratorStep )
    {
        timesVector.push_back( epoch - initialEpoch );
        resultDifferenceVector.push_back(
                earthTimeScaleConverter->getTimeDifference(
                        body_centered_coordinate_time_scale, local_proper_time_scale, epoch, stationName ) );
    }
    Eigen::VectorXd timesEigen( static_cast< int >( timesVector.size( ) ) );
    Eigen::VectorXd resultsEigen( static_cast< int >( resultDifferenceVector.size( ) ) );
    for( int i = 0; i < static_cast< int >( timesVector.size( ) ); ++i )
    {
        timesEigen( i ) = timesVector[ i ];
        resultsEigen( i ) = resultDifferenceVector[ i ];
    }
    const std::vector< double > polynomialPowers{ 0.0, 1.0 };
    Eigen::VectorXd trendFit = linear_algebra::getLeastSquaresPolynomialFit( timesEigen, resultsEigen, polynomialPowers );
    const double measuredRate = trendFit( 1 );
    std::cout<< "measuredRate " << measuredRate << std::endl;

    // For an ideal clock on the geoid, d(proper - TCG)/dt ≈ -LG (sign is negative for TCG - proper)
    BOOST_CHECK_CLOSE_FRACTION( measuredRate, -physical_constants::LG_TIME_RATE_TERM, 1.0E-6 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat

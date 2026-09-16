/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/matrixTextFileReader.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createBodyShapeModel.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::mathematical_constants;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;

BOOST_AUTO_TEST_SUITE( test_position_angle_and_separation_observation_model )

BOOST_AUTO_TEST_CASE( testPositionAngleConventionAndSingularities )
{
    const Eigen::Vector3d northPole = Eigen::Vector3d::UnitZ( );
    const Eigen::Vector3d xDirection = Eigen::Vector3d::UnitX( );
    const Eigen::Vector3d yDirection = Eigen::Vector3d::UnitY( );
    const Eigen::Vector3d negativeXDirection = -xDirection;
    const Eigen::Vector3d negativeYDirection = -yDirection;

    Eigen::Vector2d observation = calculatePositionAngleAndSeparation( xDirection, yDirection, northPole );
    BOOST_CHECK_SMALL( observation( 0 ) - mathematical_constants::PI / 2.0, 10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL( observation( 1 ) - mathematical_constants::PI / 2.0, 10.0 * std::numeric_limits< double >::epsilon( ) );

    observation = calculatePositionAngleAndSeparation( xDirection, negativeYDirection, northPole );
    BOOST_CHECK_SMALL( observation( 0 ) + mathematical_constants::PI / 2.0, 10.0 * std::numeric_limits< double >::epsilon( ) );

    BOOST_CHECK_THROW( calculatePositionAngleAndSeparation( northPole, xDirection, northPole ), std::runtime_error );
    BOOST_CHECK_THROW( calculatePositionAngleAndSeparation( xDirection, xDirection, northPole ), std::runtime_error );
    BOOST_CHECK_THROW( calculatePositionAngleAndSeparation( xDirection, negativeXDirection, northPole ), std::runtime_error );

    // Separation remains well-defined when the first line of sight is at the reference pole.
    observation = calculatePositionAngleAndSeparation( northPole, xDirection, northPole, false );
    BOOST_CHECK_SMALL( observation( 1 ) - mathematical_constants::PI / 2.0, 10.0 * std::numeric_limits< double >::epsilon( ) );
}

BOOST_AUTO_TEST_CASE( testPositionAngleAndSeparationObservationModel )
{
    spice_interface::loadStandardSpiceKernels( { paths::getSpiceKernelPath( ) + "/de430_mar097_small.bsp" } );

    // Define bodies
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Phobos" );

    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings.addSettings( "Phobos" );
    bodySettings.at( "Phobos" )->ephemerisSettings = getDefaultEphemerisSettings( "Phobos" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Define link ends
    LinkDefinition linkEnds;
    linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "" );
    linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "" );
    linkEnds[ transmitter2 ] = std::make_pair< std::string, std::string >( "Phobos", "" );

    // Create light-time correction settings
    std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back(
            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( lightTimePerturbingBodies ) );

    // Create observation settings for all three models
    std::shared_ptr< ObservationModelSettings > positionAngleSettings =
            std::make_shared< PositionAngleObservationModelSettings >( linkEnds, lightTimeCorrectionSettings );
    std::shared_ptr< ObservationModelSettings > separationSettings =
            std::make_shared< SeparationObservationModelSettings >( linkEnds, lightTimeCorrectionSettings );
    std::shared_ptr< ObservationModelSettings > positionAngleAndSeparationSettings =
            std::make_shared< PositionAngleAndSeparationObservationModelSettings >( linkEnds, lightTimeCorrectionSettings );

    // Create observation models
    std::shared_ptr< ObservationModel< 1 > > positionAngleModel =
            ObservationModelCreator< 1, double, double >::createObservationModel( positionAngleSettings, bodies );
    std::shared_ptr< ObservationModel< 1 > > separationModel =
            ObservationModelCreator< 1, double, double >::createObservationModel( separationSettings, bodies );
    std::shared_ptr< ObservationModel< 2 > > positionAngleAndSeparationModel =
            ObservationModelCreator< 2, double, double >::createObservationModel( positionAngleAndSeparationSettings, bodies );

    // Test at several epochs
    double receiverObservationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;

    std::vector< double > positionAngleLinkEndTimes;
    std::vector< double > separationLinkEndTimes;
    std::vector< double > positionAngleAndSeparationLinkEndTimes;
    std::vector< Eigen::Vector6d > positionAngleLinkEndStates;
    std::vector< Eigen::Vector6d > separationLinkEndStates;
    std::vector< Eigen::Vector6d > positionAngleAndSeparationLinkEndStates;

    Eigen::VectorXd positionAngleObservation = positionAngleModel->computeObservationsWithLinkEndData(
            receiverObservationTime, receiver, positionAngleLinkEndTimes, positionAngleLinkEndStates );
    Eigen::VectorXd separationObservation = separationModel->computeObservationsWithLinkEndData(
            receiverObservationTime, receiver, separationLinkEndTimes, separationLinkEndStates );
    Eigen::VectorXd positionAngleAndSeparationObservation = positionAngleAndSeparationModel->computeObservationsWithLinkEndData(
            receiverObservationTime, receiver, positionAngleAndSeparationLinkEndTimes, positionAngleAndSeparationLinkEndStates );

    // Verify combined model matches individual models
    BOOST_CHECK_CLOSE_FRACTION(
            positionAngleObservation( 0 ), positionAngleAndSeparationObservation( 0 ), std::numeric_limits< double >::epsilon( ) * 10.0 );
    BOOST_CHECK_CLOSE_FRACTION(
            separationObservation( 0 ), positionAngleAndSeparationObservation( 1 ), std::numeric_limits< double >::epsilon( ) * 10.0 );

    // Verify link end states match
    for( int i = 0; i < 3; i++ )
    {
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                positionAngleLinkEndStates[ i ], positionAngleAndSeparationLinkEndStates[ i ], std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                separationLinkEndStates[ i ], positionAngleAndSeparationLinkEndStates[ i ], std::numeric_limits< double >::epsilon( ) );
    }

    // ===== Cross-consistency check with RA/DEC (relative angular position) model =====
    // Compute RA/DEC from the relative angular position model's link end states,
    // then convert to position angle and separation using spherical trig.
    std::shared_ptr< ObservationModelSettings > relativeAngularPositionSettings =
            std::make_shared< ObservationModelSettings >( relative_angular_position, linkEnds, lightTimeCorrectionSettings );
    std::shared_ptr< ObservationModel< 2 > > relativeAngularPositionModel =
            ObservationModelCreator< 2, double, double >::createObservationModel( relativeAngularPositionSettings, bodies );

    std::vector< double > relativeAngularPositionLinkEndTimes;
    std::vector< Eigen::Vector6d > relativeAngularPositionLinkEndStates;
    relativeAngularPositionModel->computeObservationsWithLinkEndData(
            receiverObservationTime, receiver, relativeAngularPositionLinkEndTimes, relativeAngularPositionLinkEndStates );

    // relativeAngularPositionLinkEndStates: [0]=Mars (transmitter), [1]=Phobos (transmitter2), [2]=Earth (receiver)
    const Eigen::Matrix3d globalToJ2000 = reference_frames::getECLIPJ2000toJ2000TransformationMatrix( );
    auto computeRightAscensionAndDeclination = [ globalToJ2000 ]( const Eigen::Vector6d& transmitterState,
                                                                  const Eigen::Vector6d& receiverState ) {
        const Eigen::Vector3d relativePosition = globalToJ2000 * ( transmitterState - receiverState ).segment( 0, 3 );
        double rightAscension = 2.0 *
                std::atan( relativePosition( 1 ) /
                           ( std::sqrt( relativePosition( 0 ) * relativePosition( 0 ) + relativePosition( 1 ) * relativePosition( 1 ) ) +
                             relativePosition( 0 ) ) );
        double declination = mathematical_constants::PI / 2.0 - std::acos( relativePosition( 2 ) / relativePosition.norm( ) );
        return std::make_pair( rightAscension, declination );
    };

    auto [ marsRightAscension, marsDeclination ] =
            computeRightAscensionAndDeclination( relativeAngularPositionLinkEndStates[ 0 ], relativeAngularPositionLinkEndStates[ 2 ] );
    auto [ phobosRightAscension, phobosDeclination ] =
            computeRightAscensionAndDeclination( relativeAngularPositionLinkEndStates[ 1 ], relativeAngularPositionLinkEndStates[ 2 ] );

    // Convert RA/DEC to position angle and separation
    double rightAscensionDifference = phobosRightAscension - marsRightAscension;
    double computedSeparationDistance =
            std::acos( std::sin( marsDeclination ) * std::sin( phobosDeclination ) +
                       std::cos( marsDeclination ) * std::cos( phobosDeclination ) * std::cos( rightAscensionDifference ) );
    double computedPositionAngle =
            std::atan2( std::cos( phobosDeclination ) * std::sin( rightAscensionDifference ),
                        std::sin( phobosDeclination ) * std::cos( marsDeclination ) -
                                std::cos( phobosDeclination ) * std::sin( marsDeclination ) * std::cos( rightAscensionDifference ) );

    // Compare against dedicated models
    BOOST_CHECK_SMALL( computedPositionAngle - positionAngleObservation( 0 ), 1.0e-11 );
    BOOST_CHECK_SMALL( computedSeparationDistance - separationObservation( 0 ), 1.0e-11 );
    BOOST_CHECK_SMALL( computedPositionAngle - positionAngleAndSeparationObservation( 0 ), 1.0e-11 );
    BOOST_CHECK_SMALL( computedSeparationDistance - positionAngleAndSeparationObservation( 1 ), 1.0e-11 );

    // Test error: wrong reference link end
    BOOST_CHECK_THROW( positionAngleModel->computeObservationsWithLinkEndData(
                               receiverObservationTime, transmitter, positionAngleLinkEndTimes, positionAngleLinkEndStates ),
                       std::runtime_error );
    BOOST_CHECK_THROW(
            positionAngleAndSeparationModel->computeObservationsWithLinkEndData(
                    receiverObservationTime, transmitter, positionAngleAndSeparationLinkEndTimes, positionAngleAndSeparationLinkEndStates ),
            std::runtime_error );
}

BOOST_AUTO_TEST_CASE( testPlutoCharonHstPrefitResiduals )
{
    const std::string validationDataDirectory = paths::getTudatTestDataPath( ) + "position_angle_and_separation/";
    spice_interface::loadStandardSpiceKernels( { validationDataDirectory + "plu060_pm0001_subset.bsp" } );

    BodyListSettings bodySettings( "SSB", "J2000" );
    const std::vector< std::string > bodiesToCreate = { "Earth", "Pluto", "Charon" };
    for( const std::string& bodyName : bodiesToCreate )
    {
        bodySettings.addSettings( bodyName );
        bodySettings.at( bodyName )->ephemerisSettings = directSpiceEphemerisSettings( "SSB", "J2000", bodyName );
    }
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    LinkDefinition linkEnds;
    linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "" );
    linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Pluto", "" );
    linkEnds[ transmitter2 ] = std::make_pair< std::string, std::string >( "Charon", "" );

    const std::shared_ptr< ObservationModel< 2 > > positionAngleAndSeparationModel =
            ObservationModelCreator< 2, double, double >::createObservationModel(
                    std::make_shared< PositionAngleAndSeparationObservationModelSettings >( linkEnds ), bodies );

    // Tholen and Buie (1997), as distributed in the IMCCE Natural Satellites Data Base.
    // Columns used here are UTC Julian date at exposure midpoint, separation [arcsec], and
    // position angle [deg] measured from J2000 north through east.  Every epoch matches its
    // corresponding MAST HST exposure midpoint to within 6 ms.  Earth's centre is used as
    // receiver because the archived HST states change either observable by less than 0.002 mas.
    const Eigen::MatrixXd observations = input_output::readMatrixFromFile( validationDataDirectory + "pm0001.txt", " \t" );
    BOOST_REQUIRE_EQUAL( observations.rows( ), 60 );
    BOOST_REQUIRE( observations.cols( ) >= 3 );

    double squaredSeparationResidualSum = 0.0;
    double squaredTransversePositionAngleResidualSum = 0.0;
    double separationResidualSum = 0.0;
    double transversePositionAngleResidualSum = 0.0;
    for( Eigen::Index observationIndex = 0; observationIndex < observations.rows( ); observationIndex++ )
    {
        std::ostringstream utcJulianDate;
        utcJulianDate << "JD " << std::setprecision( 16 ) << observations( observationIndex, 0 ) << " UTC";
        const double receiverObservationTime = spice_interface::convertDateStringToEphemerisTime( utcJulianDate.str( ) );

        std::vector< double > linkEndTimes;
        std::vector< Eigen::Vector6d > linkEndStates;
        const Eigen::Vector2d computedObservation = positionAngleAndSeparationModel->computeObservationsWithLinkEndData(
                receiverObservationTime, receiver, linkEndTimes, linkEndStates );

        const double observedPositionAngle = unit_conversions::convertDegreesToRadians( observations( observationIndex, 2 ) );
        const double observedSeparation = unit_conversions::convertDegreesToRadians( observations( observationIndex, 1 ) / 3600.0 );
        const double wrappedPositionAngleResidual = std::atan2( std::sin( observedPositionAngle - computedObservation( 0 ) ),
                                                                std::cos( observedPositionAngle - computedObservation( 0 ) ) );
        const double separationResidual = observedSeparation - computedObservation( 1 );
        const double transversePositionAngleResidual = observedSeparation * wrappedPositionAngleResidual;

        squaredSeparationResidualSum += separationResidual * separationResidual;
        squaredTransversePositionAngleResidualSum += transversePositionAngleResidual * transversePositionAngleResidual;
        separationResidualSum += separationResidual;
        transversePositionAngleResidualSum += transversePositionAngleResidual;
    }

    const double radiansToMilliarcseconds = 180.0 / mathematical_constants::PI * 3600.0 * 1000.0;
    const double numberOfObservations = static_cast< double >( observations.rows( ) );
    const double separationResidualMean = separationResidualSum / numberOfObservations * radiansToMilliarcseconds;
    const double transversePositionAngleResidualMean = transversePositionAngleResidualSum / numberOfObservations * radiansToMilliarcseconds;
    const double separationResidualRms = std::sqrt( squaredSeparationResidualSum / numberOfObservations ) * radiansToMilliarcseconds;
    const double transversePositionAngleResidualRms =
            std::sqrt( squaredTransversePositionAngleResidualSum / numberOfObservations ) * radiansToMilliarcseconds;

    BOOST_TEST_MESSAGE( "pm0001 separation mean [mas]: " << std::setprecision( 15 ) << separationResidualMean );
    BOOST_TEST_MESSAGE( "pm0001 separation RMS [mas]: " << separationResidualRms );
    BOOST_TEST_MESSAGE( "pm0001 transverse PA mean [mas]: " << transversePositionAngleResidualMean );
    BOOST_TEST_MESSAGE( "pm0001 transverse PA RMS [mas]: " << transversePositionAngleResidualRms );

    // PLU060 gives milliarcsecond-level pre-fit residuals comparable to the 1997 fitted-orbit
    // residuals (2.93 mas radial and 2.70 mas transverse) stored in the source data file.
    BOOST_CHECK_SMALL( separationResidualMean - 0.536985735932, 1.0e-6 );
    BOOST_CHECK_SMALL( transversePositionAngleResidualMean - 2.048031367849, 1.0e-6 );
    BOOST_CHECK_SMALL( separationResidualRms - 4.2227973382, 1.0e-6 );
    BOOST_CHECK_SMALL( transversePositionAngleResidualRms - 4.4095016757, 1.0e-6 );
}

BOOST_AUTO_TEST_CASE( testMarsSatelliteTrueOfDatePrefitResiduals )
{
    const std::string validationDataDirectory = paths::getTudatTestDataPath( ) + "position_angle_and_separation/";
    spice_interface::loadStandardSpiceKernels( { validationDataDirectory + "mar099_de441_mm0012_subset.bsp" } );

    BodyListSettings bodySettings( "SSB", "J2000" );
    const std::vector< std::string > bodiesToCreate = { "Earth", "Mars", "Phobos", "Deimos" };
    for( const std::string& bodyName : bodiesToCreate )
    {
        bodySettings.addSettings( bodyName );
        bodySettings.at( bodyName )->ephemerisSettings = directSpiceEphemerisSettings( "SSB", "J2000", bodyName );
    }
    bodySettings.at( "Earth" )->shapeModelSettings = oblateSphericalBodyShapeSettings( 6378137.0, 1.0 / 298.257223563 );
    bodySettings.at( "Earth" )->rotationModelSettings = gcrsToItrsRotationModelSettings( basic_astrodynamics::iau_2006, "J2000" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Coordinates published by the Isaac Newton Group for the Jacobus Kapteyn Telescope.
    const Eigen::Vector3d jacobusKapteynTelescopeGeodeticPosition(
            2364.0,
            unit_conversions::convertDegreesToRadians( 28.0 + 45.0 / 60.0 + 40.1 / 3600.0 ),
            unit_conversions::convertDegreesToRadians( -( 17.0 + 52.0 / 60.0 + 41.2 / 3600.0 ) ) );
    createGroundStation( bodies.at( "Earth" ), "JKT", jacobusKapteynTelescopeGeodeticPosition, coordinate_conversions::geodetic_position );

    const std::map< int, std::pair< std::string, std::string > > observedBodyPairs = { { 1, { "Mars", "Phobos" } },
                                                                                       { 2, { "Mars", "Deimos" } },
                                                                                       { 3, { "Phobos", "Deimos" } } };
    std::map< int, std::shared_ptr< ObservationModel< 2 > > > positionAngleAndSeparationModels;
    for( const auto& observedBodyPair : observedBodyPairs )
    {
        LinkDefinition linkEnds;
        linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "JKT" );
        linkEnds[ transmitter ] = std::make_pair( observedBodyPair.second.first, std::string( "" ) );
        linkEnds[ transmitter2 ] = std::make_pair( observedBodyPair.second.second, std::string( "" ) );
        positionAngleAndSeparationModels[ observedBodyPair.first ] = ObservationModelCreator< 2, double, double >::createObservationModel(
                std::make_shared< PositionAngleAndSeparationObservationModelSettings >( linkEnds ), bodies );
    }

    const std::shared_ptr< ObservationAncillarySimulationSettings > j2000AncillarySettings =
            getPositionAngleAncillarySettings( j2000_position_angle_reference_frame );
    const std::shared_ptr< ObservationAncillarySimulationSettings > trueOfDateAncillarySettings =
            getPositionAngleAncillarySettings( true_of_date_iau_1976_1980_position_angle_reference_frame );
    const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeScaleConverter =
            earth_orientation::createDefaultTimeConverter( );

    // Jones, Sinclair, and Williams (1989), as distributed in the IMCCE Natural
    // Satellites Data Base.  Columns 14 and 15 contain separation [arcsec] and
    // position angle [deg] referred to the true equator and equinox of date.
    const Eigen::MatrixXd observations = input_output::readMatrixFromFile( validationDataDirectory + "mm0012.txt", " \t" );
    BOOST_REQUIRE_EQUAL( observations.rows( ), 166 );
    BOOST_REQUIRE_EQUAL( observations.cols( ), 18 );

    double squaredSeparationResidualSum = 0.0;
    double squaredJ2000TransversePositionAngleResidualSum = 0.0;
    double squaredTrueOfDateTransversePositionAngleResidualSum = 0.0;
    double separationResidualSum = 0.0;
    double j2000TransversePositionAngleResidualSum = 0.0;
    double trueOfDateTransversePositionAngleResidualSum = 0.0;
    double maximumSeparationModelDifference = 0.0;
    for( Eigen::Index observationIndex = 0; observationIndex < observations.rows( ); observationIndex++ )
    {
        BOOST_REQUIRE_EQUAL( static_cast< int >( observations( observationIndex, 0 ) ), 1 );
        BOOST_REQUIRE_EQUAL( static_cast< int >( observations( observationIndex, 9 ) ), 0 );
        BOOST_REQUIRE_EQUAL( static_cast< int >( observations( observationIndex, 10 ) ), 1 );
        BOOST_REQUIRE_EQUAL( static_cast< int >( observations( observationIndex, 11 ) ), 1 );
        BOOST_REQUIRE_EQUAL( static_cast< int >( observations( observationIndex, 12 ) ), 1 );

        const int year = static_cast< int >( observations( observationIndex, 2 ) );
        const int month = static_cast< int >( observations( observationIndex, 3 ) );
        const double decimalDay = observations( observationIndex, 4 );
        const int day = static_cast< int >( std::floor( decimalDay ) );
        const double secondsOfDay = ( decimalDay - static_cast< double >( day ) ) * physical_constants::JULIAN_DAY;
        const int hour = static_cast< int >( secondsOfDay / 3600.0 );
        const int minute = static_cast< int >( ( secondsOfDay - 3600.0 * hour ) / 60.0 );
        const long double second = secondsOfDay - 3600.0 * hour - 60.0 * minute;
        const basic_astrodynamics::DateTime utcDateTime( year, month, day, hour, minute, second );
        const double receiverObservationTime = timeScaleConverter->getCurrentTime(
                basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, utcDateTime.epoch< double >( ) );

        const int observedBodyPairCode = static_cast< int >( observations( observationIndex, 8 ) );
        BOOST_REQUIRE( positionAngleAndSeparationModels.count( observedBodyPairCode ) == 1 );
        std::vector< double > linkEndTimes;
        std::vector< Eigen::Vector6d > linkEndStates;
        const Eigen::Vector2d j2000ComputedObservation =
                positionAngleAndSeparationModels.at( observedBodyPairCode )
                        ->computeObservationsWithLinkEndData(
                                receiverObservationTime, receiver, linkEndTimes, linkEndStates, j2000AncillarySettings );
        const Eigen::Vector2d trueOfDateComputedObservation =
                positionAngleAndSeparationModels.at( observedBodyPairCode )
                        ->computeObservationsWithLinkEndData(
                                receiverObservationTime, receiver, linkEndTimes, linkEndStates, trueOfDateAncillarySettings );

        const double observedSeparation = unit_conversions::convertDegreesToRadians( observations( observationIndex, 13 ) / 3600.0 );
        const double observedPositionAngle = unit_conversions::convertDegreesToRadians( observations( observationIndex, 14 ) );
        const double separationResidual = observedSeparation - trueOfDateComputedObservation( 1 );
        const double j2000PositionAngleResidual = std::atan2( std::sin( observedPositionAngle - j2000ComputedObservation( 0 ) ),
                                                              std::cos( observedPositionAngle - j2000ComputedObservation( 0 ) ) );
        const double trueOfDatePositionAngleResidual = std::atan2( std::sin( observedPositionAngle - trueOfDateComputedObservation( 0 ) ),
                                                                   std::cos( observedPositionAngle - trueOfDateComputedObservation( 0 ) ) );

        squaredSeparationResidualSum += separationResidual * separationResidual;
        squaredJ2000TransversePositionAngleResidualSum += std::pow( observedSeparation * j2000PositionAngleResidual, 2 );
        squaredTrueOfDateTransversePositionAngleResidualSum += std::pow( observedSeparation * trueOfDatePositionAngleResidual, 2 );
        separationResidualSum += separationResidual;
        j2000TransversePositionAngleResidualSum += observedSeparation * j2000PositionAngleResidual;
        trueOfDateTransversePositionAngleResidualSum += observedSeparation * trueOfDatePositionAngleResidual;
        maximumSeparationModelDifference = std::max( maximumSeparationModelDifference,
                                                     std::abs( j2000ComputedObservation( 1 ) - trueOfDateComputedObservation( 1 ) ) );
    }

    const double radiansToArcseconds = 180.0 / mathematical_constants::PI * 3600.0;
    const double numberOfObservations = static_cast< double >( observations.rows( ) );
    const double separationResidualRms = std::sqrt( squaredSeparationResidualSum / numberOfObservations ) * radiansToArcseconds;
    const double j2000TransversePositionAngleResidualRms =
            std::sqrt( squaredJ2000TransversePositionAngleResidualSum / numberOfObservations ) * radiansToArcseconds;
    const double trueOfDateTransversePositionAngleResidualRms =
            std::sqrt( squaredTrueOfDateTransversePositionAngleResidualSum / numberOfObservations ) * radiansToArcseconds;
    const double separationResidualMean = separationResidualSum / numberOfObservations * radiansToArcseconds;
    const double j2000TransversePositionAngleResidualMean =
            j2000TransversePositionAngleResidualSum / numberOfObservations * radiansToArcseconds;
    const double trueOfDateTransversePositionAngleResidualMean =
            trueOfDateTransversePositionAngleResidualSum / numberOfObservations * radiansToArcseconds;

    BOOST_TEST_MESSAGE( "mm0012 separation mean [arcsec]: " << std::setprecision( 15 ) << separationResidualMean );
    BOOST_TEST_MESSAGE( "mm0012 separation RMS [arcsec]: " << separationResidualRms );
    BOOST_TEST_MESSAGE( "mm0012 J2000 transverse PA mean [arcsec]: " << j2000TransversePositionAngleResidualMean );
    BOOST_TEST_MESSAGE( "mm0012 J2000 transverse PA RMS [arcsec]: " << j2000TransversePositionAngleResidualRms );
    BOOST_TEST_MESSAGE( "mm0012 true-of-date transverse PA mean [arcsec]: " << trueOfDateTransversePositionAngleResidualMean );
    BOOST_TEST_MESSAGE( "mm0012 true-of-date transverse PA RMS [arcsec]: " << trueOfDateTransversePositionAngleResidualRms );
    BOOST_CHECK_SMALL( maximumSeparationModelDifference, 5.0 * std::numeric_limits< double >::epsilon( ) );

    // The frame correction is only about 6 mas in the transverse direction and is
    // smaller than the 0.2-arcsec measurement noise.  It therefore need not reduce
    // the unweighted RMS.  These values were independently reproduced with ERFA's
    // pnm80 matrix, SPICE light-time iteration, and the compact JPL ephemeris.
    BOOST_CHECK_SMALL( separationResidualMean + 0.020221890154, 1.0e-6 );
    BOOST_CHECK_SMALL( j2000TransversePositionAngleResidualMean - 0.041777593644, 1.0e-6 );
    BOOST_CHECK_SMALL( trueOfDateTransversePositionAngleResidualMean - 0.047758030487, 1.0e-6 );
    BOOST_CHECK_SMALL( separationResidualRms - 0.183369789212, 1.0e-6 );
    BOOST_CHECK_SMALL( j2000TransversePositionAngleResidualRms - 0.150921157572, 1.0e-6 );
    BOOST_CHECK_SMALL( trueOfDateTransversePositionAngleResidualRms - 0.152725650254, 1.0e-6 );
}

BOOST_AUTO_TEST_CASE( testSaturnSatelliteApparentDirectionPrefitResiduals )
{
    const std::string validationDataDirectory = paths::getTudatTestDataPath( ) + "position_angle_and_separation/";
    spice_interface::loadStandardSpiceKernels( { validationDataDirectory + "sat441_de441_qiao1999_subset.bsp" } );

    BodyListSettings bodySettings( "SSB", "J2000" );
    const std::vector< std::string > bodiesToCreate = { "Earth", "Enceladus", "Tethys", "Dione", "Rhea", "Titan" };
    for( const std::string& bodyName : bodiesToCreate )
    {
        bodySettings.addSettings( bodyName );
        bodySettings.at( bodyName )->ephemerisSettings = directSpiceEphemerisSettings( "SSB", "J2000", bodyName );
    }
    bodySettings.at( "Earth" )->shapeModelSettings = oblateSphericalBodyShapeSettings( 6378137.0, 1.0 / 298.257223563 );
    bodySettings.at( "Earth" )->rotationModelSettings = gcrsToItrsRotationModelSettings( basic_astrodynamics::iau_2006, "J2000" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Qiao et al. (1999) give the Sheshan observing-site coordinates. Longitude is positive east.
    const Eigen::Vector3d sheshanGeodeticPosition( 97.0,
                                                   unit_conversions::convertDegreesToRadians( 31.0 + 5.0 / 60.0 + 46.1 / 3600.0 ),
                                                   unit_conversions::convertDegreesToRadians( 121.0 + 11.0 / 60.0 + 3.3 / 3600.0 ) );
    createGroundStation( bodies.at( "Earth" ), "Sheshan", sheshanGeodeticPosition, coordinate_conversions::geodetic_position );

    const std::map< int, std::string > observedSatellites = { { 602, "Enceladus" }, { 603, "Tethys" }, { 604, "Dione" }, { 605, "Rhea" } };
    std::map< int, std::shared_ptr< ObservationModel< 2 > > > positionAngleAndSeparationModels;
    for( const auto& observedSatellite : observedSatellites )
    {
        LinkDefinition linkEnds;
        linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "Sheshan" );
        linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Titan", "" );
        linkEnds[ transmitter2 ] = std::make_pair( observedSatellite.second, std::string( "" ) );
        positionAngleAndSeparationModels[ observedSatellite.first ] = ObservationModelCreator< 2, double, double >::createObservationModel(
                std::make_shared< PositionAngleAndSeparationObservationModelSettings >( linkEnds ), bodies );
    }

    const std::shared_ptr< ObservationAncillarySimulationSettings > astrometricDirectionAncillarySettings =
            getPositionAngleAncillarySettings(
                    true_of_date_iau_1976_1980_position_angle_reference_frame, TUDAT_NAN, astrometric_position_angle_direction );
    const std::shared_ptr< ObservationAncillarySimulationSettings > aberratedDirectionAncillarySettings = getPositionAngleAncillarySettings(
            true_of_date_iau_1976_1980_position_angle_reference_frame, TUDAT_NAN, aberrated_position_angle_direction );

    // Qiao et al. (1999), Table 4. The paper labels these values apparent and topocentric,
    // with differential refraction removed but stellar aberration and parallax retained.
    // Columns are UTC Julian date, measured-satellite NAIF ID, position angle [deg], and separation [arcsec].
    const Eigen::MatrixXd observations = input_output::readMatrixFromFile( validationDataDirectory + "qiao1999_table4.txt", " \t" );
    BOOST_REQUIRE_EQUAL( observations.rows( ), 41 );
    BOOST_REQUIRE_EQUAL( observations.cols( ), 4 );

    double squaredAstrometricSeparationResidualSum = 0.0;
    double squaredAberratedSeparationResidualSum = 0.0;
    double squaredAstrometricTransversePositionAngleResidualSum = 0.0;
    double squaredAberratedTransversePositionAngleResidualSum = 0.0;
    double astrometricSeparationResidualSum = 0.0;
    double aberratedSeparationResidualSum = 0.0;
    double astrometricTransversePositionAngleResidualSum = 0.0;
    double aberratedTransversePositionAngleResidualSum = 0.0;
    double maximumSeparationAberrationCorrection = 0.0;
    double maximumTransversePositionAngleAberrationCorrection = 0.0;
    for( Eigen::Index observationIndex = 0; observationIndex < observations.rows( ); observationIndex++ )
    {
        std::ostringstream utcJulianDate;
        utcJulianDate << "JD " << std::setprecision( 16 ) << observations( observationIndex, 0 ) << " UTC";
        const double receiverObservationTime = spice_interface::convertDateStringToEphemerisTime( utcJulianDate.str( ) );
        const int observedSatelliteNaifId = static_cast< int >( observations( observationIndex, 1 ) );
        BOOST_REQUIRE( positionAngleAndSeparationModels.count( observedSatelliteNaifId ) == 1 );

        std::vector< double > linkEndTimes;
        std::vector< Eigen::Vector6d > linkEndStates;
        const Eigen::Vector2d astrometricComputedObservation =
                positionAngleAndSeparationModels.at( observedSatelliteNaifId )
                        ->computeObservationsWithLinkEndData(
                                receiverObservationTime, receiver, linkEndTimes, linkEndStates, astrometricDirectionAncillarySettings );
        const Eigen::Vector2d aberratedComputedObservation =
                positionAngleAndSeparationModels.at( observedSatelliteNaifId )
                        ->computeObservationsWithLinkEndData(
                                receiverObservationTime, receiver, linkEndTimes, linkEndStates, aberratedDirectionAncillarySettings );

        const double observedPositionAngle = unit_conversions::convertDegreesToRadians( observations( observationIndex, 2 ) );
        const double observedSeparation = unit_conversions::convertDegreesToRadians( observations( observationIndex, 3 ) / 3600.0 );
        const double astrometricPositionAngleResidual =
                std::atan2( std::sin( observedPositionAngle - astrometricComputedObservation( 0 ) ),
                            std::cos( observedPositionAngle - astrometricComputedObservation( 0 ) ) );
        const double aberratedPositionAngleResidual = std::atan2( std::sin( observedPositionAngle - aberratedComputedObservation( 0 ) ),
                                                                  std::cos( observedPositionAngle - aberratedComputedObservation( 0 ) ) );

        squaredAstrometricSeparationResidualSum += std::pow( observedSeparation - astrometricComputedObservation( 1 ), 2 );
        squaredAberratedSeparationResidualSum += std::pow( observedSeparation - aberratedComputedObservation( 1 ), 2 );
        squaredAstrometricTransversePositionAngleResidualSum += std::pow( observedSeparation * astrometricPositionAngleResidual, 2 );
        squaredAberratedTransversePositionAngleResidualSum += std::pow( observedSeparation * aberratedPositionAngleResidual, 2 );
        astrometricSeparationResidualSum += observedSeparation - astrometricComputedObservation( 1 );
        aberratedSeparationResidualSum += observedSeparation - aberratedComputedObservation( 1 );
        astrometricTransversePositionAngleResidualSum += observedSeparation * astrometricPositionAngleResidual;
        aberratedTransversePositionAngleResidualSum += observedSeparation * aberratedPositionAngleResidual;
        maximumSeparationAberrationCorrection =
                std::max( maximumSeparationAberrationCorrection,
                          std::abs( aberratedComputedObservation( 1 ) - astrometricComputedObservation( 1 ) ) );
        const double positionAngleAberrationCorrection =
                std::atan2( std::sin( aberratedComputedObservation( 0 ) - astrometricComputedObservation( 0 ) ),
                            std::cos( aberratedComputedObservation( 0 ) - astrometricComputedObservation( 0 ) ) );
        maximumTransversePositionAngleAberrationCorrection = std::max( maximumTransversePositionAngleAberrationCorrection,
                                                                       std::abs( observedSeparation * positionAngleAberrationCorrection ) );
    }

    const double radiansToArcseconds = 180.0 / mathematical_constants::PI * 3600.0;
    const double numberOfObservations = static_cast< double >( observations.rows( ) );
    const double astrometricSeparationResidualRms =
            std::sqrt( squaredAstrometricSeparationResidualSum / numberOfObservations ) * radiansToArcseconds;
    const double aberratedSeparationResidualRms =
            std::sqrt( squaredAberratedSeparationResidualSum / numberOfObservations ) * radiansToArcseconds;
    const double astrometricTransversePositionAngleResidualRms =
            std::sqrt( squaredAstrometricTransversePositionAngleResidualSum / numberOfObservations ) * radiansToArcseconds;
    const double aberratedTransversePositionAngleResidualRms =
            std::sqrt( squaredAberratedTransversePositionAngleResidualSum / numberOfObservations ) * radiansToArcseconds;
    const double astrometricSeparationResidualMean = astrometricSeparationResidualSum / numberOfObservations * radiansToArcseconds;
    const double aberratedSeparationResidualMean = aberratedSeparationResidualSum / numberOfObservations * radiansToArcseconds;
    const double astrometricTransversePositionAngleResidualMean =
            astrometricTransversePositionAngleResidualSum / numberOfObservations * radiansToArcseconds;
    const double aberratedTransversePositionAngleResidualMean =
            aberratedTransversePositionAngleResidualSum / numberOfObservations * radiansToArcseconds;

    BOOST_TEST_MESSAGE( "Qiao 1999 astrometric separation mean [arcsec]: " << std::setprecision( 15 )
                                                                           << astrometricSeparationResidualMean );
    BOOST_TEST_MESSAGE( "Qiao 1999 astrometric separation RMS [arcsec]: " << astrometricSeparationResidualRms );
    BOOST_TEST_MESSAGE( "Qiao 1999 aberrated separation mean [arcsec]: " << aberratedSeparationResidualMean );
    BOOST_TEST_MESSAGE( "Qiao 1999 aberrated separation RMS [arcsec]: " << aberratedSeparationResidualRms );
    BOOST_TEST_MESSAGE( "Qiao 1999 astrometric transverse PA mean [arcsec]: " << astrometricTransversePositionAngleResidualMean );
    BOOST_TEST_MESSAGE( "Qiao 1999 astrometric transverse PA RMS [arcsec]: " << astrometricTransversePositionAngleResidualRms );
    BOOST_TEST_MESSAGE( "Qiao 1999 aberrated transverse PA mean [arcsec]: " << aberratedTransversePositionAngleResidualMean );
    BOOST_TEST_MESSAGE( "Qiao 1999 aberrated transverse PA RMS [arcsec]: " << aberratedTransversePositionAngleResidualRms );
    BOOST_TEST_MESSAGE( "Maximum separation aberration correction [arcsec]: " << maximumSeparationAberrationCorrection *
                                radiansToArcseconds );
    BOOST_TEST_MESSAGE( "Maximum transverse PA aberration correction [arcsec]: " << maximumTransversePositionAngleAberrationCorrection *
                                radiansToArcseconds );

    // A separate SpiceyPy/ERFA calculation, using direct light-time iteration and an
    // Earth-centre receiver, reproduces the correction to within 0.12 mas. The exact
    // values below include Sheshan's position and diurnal velocity through Tudat's
    // IAU-2006 GCRS/ITRS rotation. The correction need not reduce noisy pre-fit RMS.
    BOOST_CHECK_SMALL( astrometricSeparationResidualMean + 0.094469477996, 1.0e-6 );
    BOOST_CHECK_SMALL( aberratedSeparationResidualMean + 0.100783261149, 1.0e-6 );
    BOOST_CHECK_SMALL( astrometricTransversePositionAngleResidualMean - 0.031546102027, 1.0e-6 );
    BOOST_CHECK_SMALL( aberratedTransversePositionAngleResidualMean - 0.031869600441, 1.0e-6 );
    BOOST_CHECK_SMALL( astrometricSeparationResidualRms - 0.211648598253, 1.0e-6 );
    BOOST_CHECK_SMALL( aberratedSeparationResidualRms - 0.214326363119, 1.0e-6 );
    BOOST_CHECK_SMALL( astrometricTransversePositionAngleResidualRms - 0.079196664469, 1.0e-6 );
    BOOST_CHECK_SMALL( aberratedTransversePositionAngleResidualRms - 0.079340228939, 1.0e-6 );
    BOOST_CHECK_SMALL( maximumSeparationAberrationCorrection * radiansToArcseconds - 0.007612263011, 1.0e-6 );
    BOOST_CHECK_SMALL( maximumTransversePositionAngleAberrationCorrection * radiansToArcseconds - 0.000386416656, 1.0e-6 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat

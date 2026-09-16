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

//! Check north-through-east angle signs, right-angle separation, and geometries where position angle is undefined.
BOOST_AUTO_TEST_CASE( testPositionAngleConventionAndSingularities )
{
    const Eigen::Vector3d northPole = Eigen::Vector3d::UnitZ( );
    const Eigen::Vector3d xDirection = Eigen::Vector3d::UnitX( );
    const Eigen::Vector3d yDirection = Eigen::Vector3d::UnitY( );
    const Eigen::Vector3d negativeXDirection = -xDirection;
    const Eigen::Vector3d negativeYDirection = -yDirection;

    Eigen::Vector2d observation = calculatePositionAngleAndSeparation( xDirection, yDirection, northPole );
    // A source due east of the reference direction has PA +90 degrees and separation 90 degrees.
    BOOST_CHECK_SMALL( observation( 0 ) - mathematical_constants::PI / 2.0, 10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL( observation( 1 ) - mathematical_constants::PI / 2.0, 10.0 * std::numeric_limits< double >::epsilon( ) );

    observation = calculatePositionAngleAndSeparation( xDirection, negativeYDirection, northPole );
    // Reversing east to west reverses the sign of PA.
    BOOST_CHECK_SMALL( observation( 0 ) + mathematical_constants::PI / 2.0, 10.0 * std::numeric_limits< double >::epsilon( ) );

    // PA is undefined at the pole of its tangent plane and for coincident or antipodal sight lines.
    BOOST_CHECK_THROW( calculatePositionAngleAndSeparation( northPole, xDirection, northPole ), std::runtime_error );
    BOOST_CHECK_THROW( calculatePositionAngleAndSeparation( xDirection, xDirection, northPole ), std::runtime_error );
    BOOST_CHECK_THROW( calculatePositionAngleAndSeparation( xDirection, negativeXDirection, northPole ), std::runtime_error );

    // Separation remains well-defined when the first line of sight is at the reference pole.
    observation = calculatePositionAngleAndSeparation( northPole, xDirection, northPole, false );
    BOOST_CHECK_SMALL( observation( 1 ) - mathematical_constants::PI / 2.0, 10.0 * std::numeric_limits< double >::epsilon( ) );
}

//! Compare P-only, S-only, and combined models against each other and an independent RA/Dec construction.
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

    // Simulate all three observables at one common receiver epoch with identical link ends and light-time corrections.
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

    // The two components of the combined observable must equal the corresponding standalone observables.
    BOOST_CHECK_CLOSE_FRACTION(
            positionAngleObservation( 0 ), positionAngleAndSeparationObservation( 0 ), std::numeric_limits< double >::epsilon( ) * 10.0 );
    BOOST_CHECK_CLOSE_FRACTION(
            separationObservation( 0 ), positionAngleAndSeparationObservation( 1 ), std::numeric_limits< double >::epsilon( ) * 10.0 );

    // Each wrapper must return the same retarded transmitter and receiver states as the combined model.
    for( int i = 0; i < 3; i++ )
    {
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                positionAngleLinkEndStates[ i ], positionAngleAndSeparationLinkEndStates[ i ], std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                separationLinkEndStates[ i ], positionAngleAndSeparationLinkEndStates[ i ], std::numeric_limits< double >::epsilon( ) );
    }

    // Independent geometric check: obtain retarded states from the relative-RA/Dec model, then use spherical trigonometry.
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

    // The independently reconstructed P/S values must agree with both dedicated and combined model outputs.
    BOOST_CHECK_SMALL( computedPositionAngle - positionAngleObservation( 0 ), 1.0e-11 );
    BOOST_CHECK_SMALL( computedSeparationDistance - separationObservation( 0 ), 1.0e-11 );
    BOOST_CHECK_SMALL( computedPositionAngle - positionAngleAndSeparationObservation( 0 ), 1.0e-11 );
    BOOST_CHECK_SMALL( computedSeparationDistance - positionAngleAndSeparationObservation( 1 ), 1.0e-11 );

    // P/S time tags are defined at reception; a transmitter-referenced epoch must be rejected.
    BOOST_CHECK_THROW( positionAngleModel->computeObservationsWithLinkEndData(
                               receiverObservationTime, transmitter, positionAngleLinkEndTimes, positionAngleLinkEndStates ),
                       std::runtime_error );
    BOOST_CHECK_THROW(
            positionAngleAndSeparationModel->computeObservationsWithLinkEndData(
                    receiverObservationTime, transmitter, positionAngleAndSeparationLinkEndTimes, positionAngleAndSeparationLinkEndStates ),
            std::runtime_error );
}

//! Validate default astrometric J2000 P/S against HST Pluto--Charon measurements without orbit fitting.
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
    // Ensure the expected 60-record file and required UTC/separation/PA columns are present.
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

    const double radiansToMilliArcseconds = 180.0 / mathematical_constants::PI * 3600.0 * 1000.0;
    const double numberOfObservations = static_cast< double >( observations.rows( ) );
    const double separationResidualMeanMilliArcseconds = separationResidualSum / numberOfObservations * radiansToMilliArcseconds;
    const double transversePositionAngleResidualMeanMilliArcseconds =
            transversePositionAngleResidualSum / numberOfObservations * radiansToMilliArcseconds;
    const double separationResidualRmsMilliArcseconds =
            std::sqrt( squaredSeparationResidualSum / numberOfObservations ) * radiansToMilliArcseconds;
    const double transversePositionAngleResidualRmsMilliArcseconds =
            std::sqrt( squaredTransversePositionAngleResidualSum / numberOfObservations ) * radiansToMilliArcseconds;

    BOOST_TEST_MESSAGE( "pm0001 separation mean [mas]: " << std::setprecision( 15 ) << separationResidualMeanMilliArcseconds );
    BOOST_TEST_MESSAGE( "pm0001 separation RMS [mas]: " << separationResidualRmsMilliArcseconds );
    BOOST_TEST_MESSAGE( "pm0001 transverse PA mean [mas]: " << transversePositionAngleResidualMeanMilliArcseconds );
    BOOST_TEST_MESSAGE( "pm0001 transverse PA RMS [mas]: " << transversePositionAngleResidualRmsMilliArcseconds );

    // PLU060 gives milliarcsecond-level pre-fit residuals comparable to the 1997
    // fitted-orbit residuals (2.93 mas radial and 2.70 mas transverse) stored in
    // the source file.  Check the independently expected scale rather than
    // freezing the output of the current implementation as a golden value.
    // Mean and RMS of radial and transverse residuals should have the published milliarcsecond scale.
    BOOST_CHECK_LT( std::abs( separationResidualMeanMilliArcseconds ), 5.0 );
    BOOST_CHECK_LT( std::abs( transversePositionAngleResidualMeanMilliArcseconds ), 5.0 );
    BOOST_CHECK_GT( separationResidualRmsMilliArcseconds, 1.0 );
    BOOST_CHECK_LT( separationResidualRmsMilliArcseconds, 1.0e1 );
    BOOST_CHECK_GT( transversePositionAngleResidualRmsMilliArcseconds, 1.0 );
    BOOST_CHECK_LT( transversePositionAngleResidualRmsMilliArcseconds, 1.0e1 );
}

//! Validate true-of-date PA on real Mars-satellite P/S data and isolate its effect from the invariant separation.
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
    // The record count and columns must match the documented NSDB file layout used below.
    BOOST_REQUIRE_EQUAL( observations.rows( ), 166 );
    BOOST_REQUIRE_EQUAL( observations.cols( ), 18 );

    double squaredSeparationResidualSum = 0.0;
    double squaredJ2000TransversePositionAngleResidualSum = 0.0;
    double squaredTrueOfDateTransversePositionAngleResidualSum = 0.0;
    double separationResidualSum = 0.0;
    double j2000TransversePositionAngleResidualSum = 0.0;
    double trueOfDateTransversePositionAngleResidualSum = 0.0;
    double maximumSeparationModelDifference = 0.0;
    double squaredTransversePositionAngleFrameCorrectionSum = 0.0;
    double maximumTransversePositionAngleFrameCorrection = 0.0;
    for( Eigen::Index observationIndex = 0; observationIndex < observations.rows( ); observationIndex++ )
    {
        // Restrict the comparison to the documented astrometric P/S format and available data flags.
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
        // Every coded Mars/satellite pair must map to a model with the correct transmitter ordering.
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
        const double positionAngleFrameCorrection =
                std::atan2( std::sin( trueOfDateComputedObservation( 0 ) - j2000ComputedObservation( 0 ) ),
                            std::cos( trueOfDateComputedObservation( 0 ) - j2000ComputedObservation( 0 ) ) );
        const double transversePositionAngleFrameCorrection = observedSeparation * positionAngleFrameCorrection;
        squaredTransversePositionAngleFrameCorrectionSum += transversePositionAngleFrameCorrection * transversePositionAngleFrameCorrection;
        maximumTransversePositionAngleFrameCorrection =
                std::max( maximumTransversePositionAngleFrameCorrection, std::abs( transversePositionAngleFrameCorrection ) );
    }

    const double radiansToArcseconds = 180.0 / mathematical_constants::PI * 3600.0;
    const double numberOfObservations = static_cast< double >( observations.rows( ) );
    const double separationResidualRmsArcseconds = std::sqrt( squaredSeparationResidualSum / numberOfObservations ) * radiansToArcseconds;
    const double j2000TransversePositionAngleResidualRmsArcseconds =
            std::sqrt( squaredJ2000TransversePositionAngleResidualSum / numberOfObservations ) * radiansToArcseconds;
    const double trueOfDateTransversePositionAngleResidualRmsArcseconds =
            std::sqrt( squaredTrueOfDateTransversePositionAngleResidualSum / numberOfObservations ) * radiansToArcseconds;
    const double separationResidualMeanArcseconds = separationResidualSum / numberOfObservations * radiansToArcseconds;
    const double j2000TransversePositionAngleResidualMeanArcseconds =
            j2000TransversePositionAngleResidualSum / numberOfObservations * radiansToArcseconds;
    const double trueOfDateTransversePositionAngleResidualMeanArcseconds =
            trueOfDateTransversePositionAngleResidualSum / numberOfObservations * radiansToArcseconds;
    const double transversePositionAngleFrameCorrectionRmsArcseconds =
            std::sqrt( squaredTransversePositionAngleFrameCorrectionSum / numberOfObservations ) * radiansToArcseconds;

    BOOST_TEST_MESSAGE( "mm0012 separation mean [arcsec]: " << std::setprecision( 15 ) << separationResidualMeanArcseconds );
    BOOST_TEST_MESSAGE( "mm0012 separation RMS [arcsec]: " << separationResidualRmsArcseconds );
    BOOST_TEST_MESSAGE( "mm0012 J2000 transverse PA mean [arcsec]: " << j2000TransversePositionAngleResidualMeanArcseconds );
    BOOST_TEST_MESSAGE( "mm0012 J2000 transverse PA RMS [arcsec]: " << j2000TransversePositionAngleResidualRmsArcseconds );
    BOOST_TEST_MESSAGE( "mm0012 true-of-date transverse PA mean [arcsec]: " << trueOfDateTransversePositionAngleResidualMeanArcseconds );
    BOOST_TEST_MESSAGE( "mm0012 true-of-date transverse PA RMS [arcsec]: " << trueOfDateTransversePositionAngleResidualRmsArcseconds );
    BOOST_TEST_MESSAGE( "mm0012 transverse PA frame correction RMS [arcsec]: " << transversePositionAngleFrameCorrectionRmsArcseconds );
    BOOST_TEST_MESSAGE( "mm0012 maximum transverse PA frame correction [arcsec]: " << maximumTransversePositionAngleFrameCorrection *
                                radiansToArcseconds );
    // Changing only the celestial reference pole must leave angular separation unchanged.
    BOOST_CHECK_SMALL( maximumSeparationModelDifference, 5.0 * std::numeric_limits< double >::epsilon( ) );

    // The frame correction is only about 6 mas in the transverse direction and is
    // smaller than the 0.2-arcsec measurement noise.  It therefore need not reduce
    // the unweighted RMS.  These values were independently reproduced with ERFA's
    // pnm80 matrix, SPICE light-time iteration, and the compact JPL ephemeris.
    // Validate pre-fit mean/RMS residual scales against the historical measurement precision.
    BOOST_CHECK_LT( std::abs( separationResidualMeanArcseconds ), 5.0e-2 );
    BOOST_CHECK_LT( std::abs( trueOfDateTransversePositionAngleResidualMeanArcseconds ), 1.0e-1 );
    BOOST_CHECK_GT( separationResidualRmsArcseconds, 1.0e-1 );
    BOOST_CHECK_LT( separationResidualRmsArcseconds, 5.0e-1 );
    BOOST_CHECK_GT( trueOfDateTransversePositionAngleResidualRmsArcseconds, 1.0e-1 );
    BOOST_CHECK_LT( trueOfDateTransversePositionAngleResidualRmsArcseconds, 5.0e-1 );
    // Separately require a nonzero, milliarcsecond-scale true-of-date PA correction.
    BOOST_CHECK_GT( transversePositionAngleFrameCorrectionRmsArcseconds, 1.0e-3 );
    BOOST_CHECK_LT( transversePositionAngleFrameCorrectionRmsArcseconds, 1.0e-2 );
    BOOST_CHECK_GT( maximumTransversePositionAngleFrameCorrection * radiansToArcseconds, 5.0e-3 );
    BOOST_CHECK_LT( maximumTransversePositionAngleFrameCorrection * radiansToArcseconds, 5.0e-2 );
}

//! Validate B1950 north on historical Nereid astrometry and resolve its effect on PA separately from residual noise.
BOOST_AUTO_TEST_CASE( testNereidB1950ReferencePolePrefitResiduals )
{
    const std::string validationDataDirectory = paths::getTudatTestDataPath( ) + "position_angle_and_separation/";
    spice_interface::loadStandardSpiceKernels( { validationDataDirectory + "nep105_nm1007_subset.bsp" } );

    BodyListSettings bodySettings( "SSB", "J2000" );
    const std::vector< std::string > bodiesToCreate = { "Earth", "Neptune", "Nereid" };
    for( const std::string& bodyName : bodiesToCreate )
    {
        bodySettings.addSettings( bodyName );
        bodySettings.at( bodyName )->ephemerisSettings = directSpiceEphemerisSettings( "SSB", "J2000", bodyName );
    }
    bodySettings.at( "Earth" )->shapeModelSettings = oblateSphericalBodyShapeSettings( 6378137.0, 1.0 / 298.257223563 );
    bodySettings.at( "Earth" )->rotationModelSettings = gcrsToItrsRotationModelSettings( basic_astrodynamics::iau_2006, "J2000" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // MPC parallax constants resolve the exact historical sites without
    // introducing a geodetic-datum interpretation.  Coordinates are Earth-fixed.
    const double terrestrialEquatorialRadius = 6378137.0;
    const std::map< int, Eigen::Vector3d > stationParallaxConstants = {
        { 568, Eigen::Vector3d( unit_conversions::convertDegreesToRadians( 204.5278 ), 0.94171, 0.33725 ) },
        { 711, Eigen::Vector3d( unit_conversions::convertDegreesToRadians( 255.9785 ), 0.86114, 0.50731 ) },
        { 809, Eigen::Vector3d( unit_conversions::convertDegreesToRadians( 289.26626 ), 0.873440, -0.486052 ) }
    };
    for( const auto& stationParallaxConstant : stationParallaxConstants )
    {
        const double longitude = stationParallaxConstant.second( 0 );
        const double distanceFromSpinAxis = terrestrialEquatorialRadius * stationParallaxConstant.second( 1 );
        const Eigen::Vector3d bodyFixedStationPosition( distanceFromSpinAxis * std::cos( longitude ),
                                                        distanceFromSpinAxis * std::sin( longitude ),
                                                        terrestrialEquatorialRadius * stationParallaxConstant.second( 2 ) );
        createGroundStation( bodies.at( "Earth" ),
                             std::to_string( stationParallaxConstant.first ),
                             bodyFixedStationPosition,
                             coordinate_conversions::cartesian_position );
    }

    std::map< int, std::shared_ptr< ObservationModel< 2 > > > positionAngleAndSeparationModels;
    for( const auto& stationParallaxConstant : stationParallaxConstants )
    {
        LinkDefinition linkEnds;
        linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth", std::to_string( stationParallaxConstant.first ) );
        linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Neptune", "" );
        linkEnds[ transmitter2 ] = std::make_pair< std::string, std::string >( "Nereid", "" );
        positionAngleAndSeparationModels[ stationParallaxConstant.first ] =
                ObservationModelCreator< 2, double, double >::createObservationModel(
                        std::make_shared< PositionAngleAndSeparationObservationModelSettings >( linkEnds ), bodies );
    }

    const std::shared_ptr< ObservationAncillarySimulationSettings > j2000AncillarySettings =
            getPositionAngleAncillarySettings( j2000_position_angle_reference_frame );
    const std::shared_ptr< ObservationAncillarySimulationSettings > b1950AncillarySettings =
            getPositionAngleAncillarySettings( b1950_position_angle_reference_frame );
    const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeScaleConverter =
            earth_orientation::createDefaultTimeConverter( );

    // Veillet (1982) and Veillet & Bois (1988), NSDB nm1007.  NSDB specifies
    // topocentric astrometric B1950 coordinates and UTC.  Atmospheric reduction
    // is undocumented and deliberately not modeled in this frame-focused test.
    const Eigen::MatrixXd observations = input_output::readMatrixFromFile( validationDataDirectory + "nm1007.txt", " \t" );
    // Check the number of usable records and the fixed-format columns read in this test.
    BOOST_REQUIRE_EQUAL( observations.rows( ), 17 );
    BOOST_REQUIRE_EQUAL( observations.cols( ), 10 );

    double separationResidualSum = 0.0;
    double squaredSeparationResidualSum = 0.0;
    double b1950TransversePositionAngleResidualSum = 0.0;
    double squaredB1950TransversePositionAngleResidualSum = 0.0;
    double squaredJ2000TransversePositionAngleResidualSum = 0.0;
    double squaredTransversePositionAngleFrameCorrectionSum = 0.0;
    double maximumTransversePositionAngleFrameCorrection = 0.0;
    double maximumSeparationFrameDifference = 0.0;
    for( Eigen::Index observationIndex = 0; observationIndex < observations.rows( ); observationIndex++ )
    {
        const int year = static_cast< int >( observations( observationIndex, 1 ) );
        const int month = static_cast< int >( observations( observationIndex, 2 ) );
        const double decimalDay = observations( observationIndex, 3 );
        const int day = static_cast< int >( std::floor( decimalDay ) );
        const double secondsOfDay = ( decimalDay - static_cast< double >( day ) ) * physical_constants::JULIAN_DAY;
        const int hour = static_cast< int >( secondsOfDay / 3600.0 );
        const int minute = static_cast< int >( ( secondsOfDay - 3600.0 * hour ) / 60.0 );
        const long double second = secondsOfDay - 3600.0 * hour - 60.0 * minute;
        const basic_astrodynamics::DateTime utcDateTime( year, month, day, hour, minute, second );
        const double receiverObservationTime = timeScaleConverter->getCurrentTime(
                basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, utcDateTime.epoch< double >( ) );

        const int observatoryCode = static_cast< int >( observations( observationIndex, 6 ) );
        // The historical observatory code must select one of the three reconstructed topocentric receivers.
        BOOST_REQUIRE( positionAngleAndSeparationModels.count( observatoryCode ) == 1 );
        std::vector< double > linkEndTimes;
        std::vector< Eigen::Vector6d > linkEndStates;
        const Eigen::Vector2d j2000ComputedObservation =
                positionAngleAndSeparationModels.at( observatoryCode )
                        ->computeObservationsWithLinkEndData(
                                receiverObservationTime, receiver, linkEndTimes, linkEndStates, j2000AncillarySettings );
        const Eigen::Vector2d b1950ComputedObservation =
                positionAngleAndSeparationModels.at( observatoryCode )
                        ->computeObservationsWithLinkEndData(
                                receiverObservationTime, receiver, linkEndTimes, linkEndStates, b1950AncillarySettings );

        const double observedX = unit_conversions::convertArcSecondsToRadians( observations( observationIndex, 4 ) );
        const double observedY = unit_conversions::convertArcSecondsToRadians( observations( observationIndex, 5 ) );
        const double observedPositionAngle = std::atan2( observedX, observedY );
        const double observedSeparation = std::hypot( observedX, observedY );
        const double b1950PositionAngleResidual = std::atan2( std::sin( observedPositionAngle - b1950ComputedObservation( 0 ) ),
                                                              std::cos( observedPositionAngle - b1950ComputedObservation( 0 ) ) );
        const double j2000PositionAngleResidual = std::atan2( std::sin( observedPositionAngle - j2000ComputedObservation( 0 ) ),
                                                              std::cos( observedPositionAngle - j2000ComputedObservation( 0 ) ) );
        const double separationResidual = observedSeparation - b1950ComputedObservation( 1 );
        const double transversePositionAngleFrameCorrection = observedSeparation *
                std::atan2( std::sin( b1950ComputedObservation( 0 ) - j2000ComputedObservation( 0 ) ),
                            std::cos( b1950ComputedObservation( 0 ) - j2000ComputedObservation( 0 ) ) );

        separationResidualSum += separationResidual;
        squaredSeparationResidualSum += separationResidual * separationResidual;
        b1950TransversePositionAngleResidualSum += observedSeparation * b1950PositionAngleResidual;
        squaredB1950TransversePositionAngleResidualSum += std::pow( observedSeparation * b1950PositionAngleResidual, 2 );
        squaredJ2000TransversePositionAngleResidualSum += std::pow( observedSeparation * j2000PositionAngleResidual, 2 );
        squaredTransversePositionAngleFrameCorrectionSum += transversePositionAngleFrameCorrection * transversePositionAngleFrameCorrection;
        maximumTransversePositionAngleFrameCorrection =
                std::max( maximumTransversePositionAngleFrameCorrection, std::abs( transversePositionAngleFrameCorrection ) );
        maximumSeparationFrameDifference =
                std::max( maximumSeparationFrameDifference, std::abs( b1950ComputedObservation( 1 ) - j2000ComputedObservation( 1 ) ) );
    }

    const double radiansToArcseconds = 180.0 / mathematical_constants::PI * 3600.0;
    const double numberOfObservations = static_cast< double >( observations.rows( ) );
    const double separationResidualMeanArcseconds = separationResidualSum / numberOfObservations * radiansToArcseconds;
    const double separationResidualRmsArcseconds = std::sqrt( squaredSeparationResidualSum / numberOfObservations ) * radiansToArcseconds;
    const double b1950TransversePositionAngleResidualMeanArcseconds =
            b1950TransversePositionAngleResidualSum / numberOfObservations * radiansToArcseconds;
    const double b1950TransversePositionAngleResidualRmsArcseconds =
            std::sqrt( squaredB1950TransversePositionAngleResidualSum / numberOfObservations ) * radiansToArcseconds;
    const double j2000TransversePositionAngleResidualRmsArcseconds =
            std::sqrt( squaredJ2000TransversePositionAngleResidualSum / numberOfObservations ) * radiansToArcseconds;
    const double transversePositionAngleFrameCorrectionRmsArcseconds =
            std::sqrt( squaredTransversePositionAngleFrameCorrectionSum / numberOfObservations ) * radiansToArcseconds;

    BOOST_TEST_MESSAGE( "nm1007 B1950 separation mean [arcsec]: " << std::setprecision( 15 ) << separationResidualMeanArcseconds );
    BOOST_TEST_MESSAGE( "nm1007 B1950 separation RMS [arcsec]: " << separationResidualRmsArcseconds );
    BOOST_TEST_MESSAGE( "nm1007 B1950 transverse PA mean [arcsec]: " << b1950TransversePositionAngleResidualMeanArcseconds );
    BOOST_TEST_MESSAGE( "nm1007 B1950 transverse PA RMS [arcsec]: " << b1950TransversePositionAngleResidualRmsArcseconds );
    BOOST_TEST_MESSAGE( "nm1007 J2000 transverse PA RMS [arcsec]: " << j2000TransversePositionAngleResidualRmsArcseconds );
    BOOST_TEST_MESSAGE( "nm1007 transverse PA frame correction RMS [arcsec]: " << transversePositionAngleFrameCorrectionRmsArcseconds );
    BOOST_TEST_MESSAGE( "nm1007 maximum transverse PA frame correction [arcsec]: " << maximumTransversePositionAngleFrameCorrection *
                                radiansToArcseconds );

    // Yuan et al. (2021) report 0.17--0.22 arcsec coordinate RMS for these
    // historical series.  The B1950 pole signal is independently much larger
    // than both the mean and the residual scale, and separation is pole-invariant.
    // Bound the B1950 pre-fit radial and transverse mean/RMS at the published historical precision.
    BOOST_CHECK_LT( std::abs( separationResidualMeanArcseconds ), 1.0e-1 );
    BOOST_CHECK_LT( std::abs( b1950TransversePositionAngleResidualMeanArcseconds ), 2.0e-1 );
    BOOST_CHECK_GT( separationResidualRmsArcseconds, 5.0e-2 );
    BOOST_CHECK_LT( separationResidualRmsArcseconds, 5.0e-1 );
    BOOST_CHECK_GT( b1950TransversePositionAngleResidualRmsArcseconds, 1.0e-1 );
    BOOST_CHECK_LT( b1950TransversePositionAngleResidualRmsArcseconds, 5.0e-1 );
    // The B1950 pole should outperform J2000 on these B1950 data and produce a resolvable PA correction.
    BOOST_CHECK_GT( j2000TransversePositionAngleResidualRmsArcseconds, b1950TransversePositionAngleResidualRmsArcseconds );
    BOOST_CHECK_GT( transversePositionAngleFrameCorrectionRmsArcseconds, 2.0e-1 );
    BOOST_CHECK_LT( transversePositionAngleFrameCorrectionRmsArcseconds, 5.0e-1 );
    BOOST_CHECK_GT( maximumTransversePositionAngleFrameCorrection * radiansToArcseconds, 5.0e-1 );
    BOOST_CHECK_LT( maximumTransversePositionAngleFrameCorrection * radiansToArcseconds, 1.0 );
    // Separation depends on the two sight lines, not on the chosen celestial north pole.
    BOOST_CHECK_SMALL( maximumSeparationFrameDifference, 5.0 * std::numeric_limits< double >::epsilon( ) );
}

//! Validate native B1950 P/S measurements of Uranian satellites and pole-invariance of separation.
BOOST_AUTO_TEST_CASE( testUranianSatelliteNativeB1950PrefitResiduals )
{
    const std::string validationDataDirectory = paths::getTudatTestDataPath( ) + "position_angle_and_separation/";
    spice_interface::loadStandardSpiceKernels( { validationDataDirectory + "ura184_de441_um0027_subset.bsp" } );

    BodyListSettings bodySettings( "SSB", "J2000" );
    const std::vector< std::string > bodiesToCreate = { "Earth", "Uranus", "Ariel", "Umbriel", "Titania", "Oberon" };
    for( const std::string& bodyName : bodiesToCreate )
    {
        bodySettings.addSettings( bodyName );
        bodySettings.at( bodyName )->ephemerisSettings = directSpiceEphemerisSettings( "SSB", "J2000", bodyName );
    }
    bodySettings.at( "Earth" )->shapeModelSettings = oblateSphericalBodyShapeSettings( 6378137.0, 1.0 / 298.257223563 );
    bodySettings.at( "Earth" )->rotationModelSettings = gcrsToItrsRotationModelSettings( basic_astrodynamics::iau_2006, "J2000" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // MPC observatory code 371: Tokyo-Okayama.  The parallax constants define
    // its Earth-fixed position without requiring an assumption about geodetic datum.
    const double terrestrialEquatorialRadius = 6378137.0;
    const double observatoryLongitude = unit_conversions::convertDegreesToRadians( 133.5965 );
    const double observatoryDistanceFromSpinAxis = terrestrialEquatorialRadius * 0.82433;
    const Eigen::Vector3d bodyFixedObservatoryPosition( observatoryDistanceFromSpinAxis * std::cos( observatoryLongitude ),
                                                        observatoryDistanceFromSpinAxis * std::sin( observatoryLongitude ),
                                                        terrestrialEquatorialRadius * 0.56431 );
    createGroundStation( bodies.at( "Earth" ), "Tokyo-Okayama", bodyFixedObservatoryPosition, coordinate_conversions::cartesian_position );

    const std::map< int, std::string > observedSatellites = { { 1, "Ariel" }, { 2, "Umbriel" }, { 3, "Titania" }, { 4, "Oberon" } };
    std::map< int, std::shared_ptr< ObservationModel< 2 > > > positionAngleAndSeparationModels;
    for( const auto& observedSatellite : observedSatellites )
    {
        LinkDefinition linkEnds;
        linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "Tokyo-Okayama" );
        linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Uranus", "" );
        linkEnds[ transmitter2 ] = std::make_pair( observedSatellite.second, std::string( "" ) );
        positionAngleAndSeparationModels[ observedSatellite.first ] = ObservationModelCreator< 2, double, double >::createObservationModel(
                std::make_shared< PositionAngleAndSeparationObservationModelSettings >( linkEnds ), bodies );
    }

    const std::shared_ptr< ObservationAncillarySimulationSettings > j2000AncillarySettings =
            getPositionAngleAncillarySettings( j2000_position_angle_reference_frame );
    const std::shared_ptr< ObservationAncillarySimulationSettings > b1950AncillarySettings =
            getPositionAngleAncillarySettings( b1950_position_angle_reference_frame );

    // Tomita and Soma (1979), NSDB um0027: native separation/position-angle
    // observations explicitly documented as topocentric, astrometric, B1950,
    // and ET.  No RA/declination-to-P/S conversion is performed here.
    const Eigen::MatrixXd observations = input_output::readMatrixFromFile( validationDataDirectory + "um0027.txt", " \t" );
    // Confirm all native P/S records and the date, body-code, separation, and PA columns are available.
    BOOST_REQUIRE_EQUAL( observations.rows( ), 156 );
    BOOST_REQUIRE_EQUAL( observations.cols( ), 10 );

    double separationResidualSum = 0.0;
    double squaredSeparationResidualSum = 0.0;
    double b1950TransversePositionAngleResidualSum = 0.0;
    double squaredB1950TransversePositionAngleResidualSum = 0.0;
    double squaredJ2000TransversePositionAngleResidualSum = 0.0;
    double squaredTransversePositionAngleFrameCorrectionSum = 0.0;
    double maximumTransversePositionAngleFrameCorrection = 0.0;
    double maximumSeparationFrameDifference = 0.0;
    for( Eigen::Index observationIndex = 0; observationIndex < observations.rows( ); observationIndex++ )
    {
        const int year = static_cast< int >( observations( observationIndex, 0 ) );
        const int month = static_cast< int >( observations( observationIndex, 1 ) );
        const int day = static_cast< int >( observations( observationIndex, 2 ) );
        const int hour = static_cast< int >( observations( observationIndex, 3 ) );
        const int minute = static_cast< int >( observations( observationIndex, 4 ) );
        const long double second = static_cast< long double >( observations( observationIndex, 5 ) );
        // NSDB specifies ET.  DateTime::epoch is the required uniform seconds
        // since J2000 representation; historical ET--TDB differences are far
        // below the precision of these photographic observations.
        const double receiverObservationTime = basic_astrodynamics::DateTime( year, month, day, hour, minute, second ).epoch< double >( );

        const int observedSatelliteCode = static_cast< int >( observations( observationIndex, 7 ) );
        // Each satellite code must select the corresponding Uranus-to-moon observation model.
        BOOST_REQUIRE( positionAngleAndSeparationModels.count( observedSatelliteCode ) == 1 );
        std::vector< double > linkEndTimes;
        std::vector< Eigen::Vector6d > linkEndStates;
        const Eigen::Vector2d j2000ComputedObservation =
                positionAngleAndSeparationModels.at( observedSatelliteCode )
                        ->computeObservationsWithLinkEndData(
                                receiverObservationTime, receiver, linkEndTimes, linkEndStates, j2000AncillarySettings );
        const Eigen::Vector2d b1950ComputedObservation =
                positionAngleAndSeparationModels.at( observedSatelliteCode )
                        ->computeObservationsWithLinkEndData(
                                receiverObservationTime, receiver, linkEndTimes, linkEndStates, b1950AncillarySettings );

        const double observedSeparation = unit_conversions::convertArcSecondsToRadians( observations( observationIndex, 8 ) );
        const double observedPositionAngle = unit_conversions::convertDegreesToRadians( observations( observationIndex, 9 ) );
        const double b1950PositionAngleResidual = std::atan2( std::sin( observedPositionAngle - b1950ComputedObservation( 0 ) ),
                                                              std::cos( observedPositionAngle - b1950ComputedObservation( 0 ) ) );
        const double j2000PositionAngleResidual = std::atan2( std::sin( observedPositionAngle - j2000ComputedObservation( 0 ) ),
                                                              std::cos( observedPositionAngle - j2000ComputedObservation( 0 ) ) );
        const double separationResidual = observedSeparation - b1950ComputedObservation( 1 );
        const double transversePositionAngleFrameCorrection = observedSeparation *
                std::atan2( std::sin( b1950ComputedObservation( 0 ) - j2000ComputedObservation( 0 ) ),
                            std::cos( b1950ComputedObservation( 0 ) - j2000ComputedObservation( 0 ) ) );

        separationResidualSum += separationResidual;
        squaredSeparationResidualSum += separationResidual * separationResidual;
        b1950TransversePositionAngleResidualSum += observedSeparation * b1950PositionAngleResidual;
        squaredB1950TransversePositionAngleResidualSum += std::pow( observedSeparation * b1950PositionAngleResidual, 2 );
        squaredJ2000TransversePositionAngleResidualSum += std::pow( observedSeparation * j2000PositionAngleResidual, 2 );
        squaredTransversePositionAngleFrameCorrectionSum += transversePositionAngleFrameCorrection * transversePositionAngleFrameCorrection;
        maximumTransversePositionAngleFrameCorrection =
                std::max( maximumTransversePositionAngleFrameCorrection, std::abs( transversePositionAngleFrameCorrection ) );
        maximumSeparationFrameDifference =
                std::max( maximumSeparationFrameDifference, std::abs( b1950ComputedObservation( 1 ) - j2000ComputedObservation( 1 ) ) );
    }

    const double radiansToArcseconds = 180.0 / mathematical_constants::PI * 3600.0;
    const double numberOfObservations = static_cast< double >( observations.rows( ) );
    const double separationResidualMeanArcseconds = separationResidualSum / numberOfObservations * radiansToArcseconds;
    const double separationResidualRmsArcseconds = std::sqrt( squaredSeparationResidualSum / numberOfObservations ) * radiansToArcseconds;
    const double b1950TransversePositionAngleResidualMeanArcseconds =
            b1950TransversePositionAngleResidualSum / numberOfObservations * radiansToArcseconds;
    const double b1950TransversePositionAngleResidualRmsArcseconds =
            std::sqrt( squaredB1950TransversePositionAngleResidualSum / numberOfObservations ) * radiansToArcseconds;
    const double j2000TransversePositionAngleResidualRmsArcseconds =
            std::sqrt( squaredJ2000TransversePositionAngleResidualSum / numberOfObservations ) * radiansToArcseconds;
    const double transversePositionAngleFrameCorrectionRmsArcseconds =
            std::sqrt( squaredTransversePositionAngleFrameCorrectionSum / numberOfObservations ) * radiansToArcseconds;

    BOOST_TEST_MESSAGE( "um0027 B1950 separation mean [arcsec]: " << std::setprecision( 15 ) << separationResidualMeanArcseconds );
    BOOST_TEST_MESSAGE( "um0027 B1950 separation RMS [arcsec]: " << separationResidualRmsArcseconds );
    BOOST_TEST_MESSAGE( "um0027 B1950 transverse PA mean [arcsec]: " << b1950TransversePositionAngleResidualMeanArcseconds );
    BOOST_TEST_MESSAGE( "um0027 B1950 transverse PA RMS [arcsec]: " << b1950TransversePositionAngleResidualRmsArcseconds );
    BOOST_TEST_MESSAGE( "um0027 J2000 transverse PA RMS [arcsec]: " << j2000TransversePositionAngleResidualRmsArcseconds );
    BOOST_TEST_MESSAGE( "um0027 transverse PA frame correction RMS [arcsec]: " << transversePositionAngleFrameCorrectionRmsArcseconds );
    BOOST_TEST_MESSAGE( "um0027 maximum transverse PA frame correction [arcsec]: " << maximumTransversePositionAngleFrameCorrection *
                                radiansToArcseconds );

    // Jacobson (2014), Table 3 reports post-fit RMS values of 0.443 arcsec in
    // separation and 0.363 arcsec transverse in position angle for retained
    // Tomita--Soma observations.  This unfiltered NSDB set has additional
    // outliers but the same sub-arcsecond-to-arcsecond scale.  The B1950
    // correction itself is at the 0.05-arcsec RMS scale.  Bounds are
    // intentionally decade-style physical checks, not frozen computed values.
    // Check that B1950 pre-fit radial and transverse means/RMS are on the historical arcsecond scale.
    BOOST_CHECK_LT( std::abs( separationResidualMeanArcseconds ), 2.0e-1 );
    BOOST_CHECK_LT( std::abs( b1950TransversePositionAngleResidualMeanArcseconds ), 2.0e-1 );
    BOOST_CHECK_GT( separationResidualRmsArcseconds, 2.0e-1 );
    BOOST_CHECK_LT( separationResidualRmsArcseconds, 2.0 );
    BOOST_CHECK_GT( b1950TransversePositionAngleResidualRmsArcseconds, 2.0e-1 );
    BOOST_CHECK_LT( b1950TransversePositionAngleResidualRmsArcseconds, 2.0 );
    // On explicitly B1950 data, J2000 should worsen PA RMS and the frame correction should be measurable.
    BOOST_CHECK_GT( j2000TransversePositionAngleResidualRmsArcseconds, b1950TransversePositionAngleResidualRmsArcseconds );
    BOOST_CHECK_GT( transversePositionAngleFrameCorrectionRmsArcseconds, 2.0e-2 );
    BOOST_CHECK_LT( transversePositionAngleFrameCorrectionRmsArcseconds, 1.0e-1 );
    BOOST_CHECK_GT( maximumTransversePositionAngleFrameCorrection * radiansToArcseconds, 5.0e-2 );
    BOOST_CHECK_LT( maximumTransversePositionAngleFrameCorrection * radiansToArcseconds, 2.0e-1 );
    // Switching north poles must not change the separation returned by the same line-of-sight geometry.
    BOOST_CHECK_SMALL( maximumSeparationFrameDifference, 5.0 * std::numeric_limits< double >::epsilon( ) );
}

//! Validate the apparent-direction option on Saturnian P/S data and resolve differential stellar-aberration effects.
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
    // Confirm the 41 records and four columns: UTC Julian date, satellite ID, PA, and separation.
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
    double squaredSeparationAberrationCorrectionSum = 0.0;
    double squaredTransversePositionAngleAberrationCorrectionSum = 0.0;
    double maximumSeparationAberrationCorrection = 0.0;
    double maximumTransversePositionAngleAberrationCorrection = 0.0;
    for( Eigen::Index observationIndex = 0; observationIndex < observations.rows( ); observationIndex++ )
    {
        std::ostringstream utcJulianDate;
        utcJulianDate << "JD " << std::setprecision( 16 ) << observations( observationIndex, 0 ) << " UTC";
        const double receiverObservationTime = spice_interface::convertDateStringToEphemerisTime( utcJulianDate.str( ) );
        const int observedSatelliteNaifId = static_cast< int >( observations( observationIndex, 1 ) );
        // The measured satellite ID must map to a Titan-reference P/S model.
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
        const double separationAberrationCorrection = aberratedComputedObservation( 1 ) - astrometricComputedObservation( 1 );
        squaredSeparationAberrationCorrectionSum += separationAberrationCorrection * separationAberrationCorrection;
        maximumSeparationAberrationCorrection =
                std::max( maximumSeparationAberrationCorrection, std::abs( separationAberrationCorrection ) );
        const double positionAngleAberrationCorrection =
                std::atan2( std::sin( aberratedComputedObservation( 0 ) - astrometricComputedObservation( 0 ) ),
                            std::cos( aberratedComputedObservation( 0 ) - astrometricComputedObservation( 0 ) ) );
        const double transversePositionAngleAberrationCorrection = observedSeparation * positionAngleAberrationCorrection;
        squaredTransversePositionAngleAberrationCorrectionSum +=
                transversePositionAngleAberrationCorrection * transversePositionAngleAberrationCorrection;
        maximumTransversePositionAngleAberrationCorrection =
                std::max( maximumTransversePositionAngleAberrationCorrection, std::abs( transversePositionAngleAberrationCorrection ) );
    }

    const double radiansToArcseconds = 180.0 / mathematical_constants::PI * 3600.0;
    const double numberOfObservations = static_cast< double >( observations.rows( ) );
    const double astrometricSeparationResidualRmsArcseconds =
            std::sqrt( squaredAstrometricSeparationResidualSum / numberOfObservations ) * radiansToArcseconds;
    const double aberratedSeparationResidualRmsArcseconds =
            std::sqrt( squaredAberratedSeparationResidualSum / numberOfObservations ) * radiansToArcseconds;
    const double astrometricTransversePositionAngleResidualRmsArcseconds =
            std::sqrt( squaredAstrometricTransversePositionAngleResidualSum / numberOfObservations ) * radiansToArcseconds;
    const double aberratedTransversePositionAngleResidualRmsArcseconds =
            std::sqrt( squaredAberratedTransversePositionAngleResidualSum / numberOfObservations ) * radiansToArcseconds;
    const double astrometricSeparationResidualMeanArcseconds =
            astrometricSeparationResidualSum / numberOfObservations * radiansToArcseconds;
    const double aberratedSeparationResidualMeanArcseconds = aberratedSeparationResidualSum / numberOfObservations * radiansToArcseconds;
    const double astrometricTransversePositionAngleResidualMeanArcseconds =
            astrometricTransversePositionAngleResidualSum / numberOfObservations * radiansToArcseconds;
    const double aberratedTransversePositionAngleResidualMeanArcseconds =
            aberratedTransversePositionAngleResidualSum / numberOfObservations * radiansToArcseconds;
    const double separationAberrationCorrectionRmsArcseconds =
            std::sqrt( squaredSeparationAberrationCorrectionSum / numberOfObservations ) * radiansToArcseconds;
    const double transversePositionAngleAberrationCorrectionRmsArcseconds =
            std::sqrt( squaredTransversePositionAngleAberrationCorrectionSum / numberOfObservations ) * radiansToArcseconds;

    BOOST_TEST_MESSAGE( "Qiao 1999 astrometric separation mean [arcsec]: " << std::setprecision( 15 )
                                                                           << astrometricSeparationResidualMeanArcseconds );
    BOOST_TEST_MESSAGE( "Qiao 1999 astrometric separation RMS [arcsec]: " << astrometricSeparationResidualRmsArcseconds );
    BOOST_TEST_MESSAGE( "Qiao 1999 aberrated separation mean [arcsec]: " << aberratedSeparationResidualMeanArcseconds );
    BOOST_TEST_MESSAGE( "Qiao 1999 aberrated separation RMS [arcsec]: " << aberratedSeparationResidualRmsArcseconds );
    BOOST_TEST_MESSAGE( "Qiao 1999 astrometric transverse PA mean [arcsec]: " << astrometricTransversePositionAngleResidualMeanArcseconds );
    BOOST_TEST_MESSAGE( "Qiao 1999 astrometric transverse PA RMS [arcsec]: " << astrometricTransversePositionAngleResidualRmsArcseconds );
    BOOST_TEST_MESSAGE( "Qiao 1999 aberrated transverse PA mean [arcsec]: " << aberratedTransversePositionAngleResidualMeanArcseconds );
    BOOST_TEST_MESSAGE( "Qiao 1999 aberrated transverse PA RMS [arcsec]: " << aberratedTransversePositionAngleResidualRmsArcseconds );
    BOOST_TEST_MESSAGE( "Separation aberration correction RMS [arcsec]: " << separationAberrationCorrectionRmsArcseconds );
    BOOST_TEST_MESSAGE( "Transverse PA aberration correction RMS [arcsec]: " << transversePositionAngleAberrationCorrectionRmsArcseconds );
    BOOST_TEST_MESSAGE( "Maximum separation aberration correction [arcsec]: " << maximumSeparationAberrationCorrection *
                                radiansToArcseconds );
    BOOST_TEST_MESSAGE( "Maximum transverse PA aberration correction [arcsec]: " << maximumTransversePositionAngleAberrationCorrection *
                                radiansToArcseconds );

    // Qiao et al. report roughly 0.08--0.20 arcsec RMS, depending on satellite
    // and coordinate.  A separate SpiceyPy calculation confirms a differential
    // aberration signal of several milliarcseconds.  Check those physical scales,
    // without making the current numerical outputs their own reference truth.
    // Compare the apparent-model pre-fit mean/RMS with the reported sub-arcsecond scatter.
    BOOST_CHECK_LT( std::abs( aberratedSeparationResidualMeanArcseconds ), 2.0e-1 );
    BOOST_CHECK_LT( std::abs( aberratedTransversePositionAngleResidualMeanArcseconds ), 1.0e-1 );
    BOOST_CHECK_GT( aberratedSeparationResidualRmsArcseconds, 1.0e-1 );
    BOOST_CHECK_LT( aberratedSeparationResidualRmsArcseconds, 5.0e-1 );
    BOOST_CHECK_GT( aberratedTransversePositionAngleResidualRmsArcseconds, 5.0e-2 );
    BOOST_CHECK_LT( aberratedTransversePositionAngleResidualRmsArcseconds, 2.0e-1 );
    // Require nonzero milliarcsecond-scale differential aberration in S and transverse PA.
    BOOST_CHECK_GT( separationAberrationCorrectionRmsArcseconds, 2.0e-3 );
    BOOST_CHECK_LT( separationAberrationCorrectionRmsArcseconds, 1.0e-2 );
    BOOST_CHECK_GT( transversePositionAngleAberrationCorrectionRmsArcseconds, 1.0e-4 );
    BOOST_CHECK_LT( transversePositionAngleAberrationCorrectionRmsArcseconds, 1.0e-3 );
    BOOST_CHECK_GT( maximumSeparationAberrationCorrection * radiansToArcseconds, 5.0e-3 );
    BOOST_CHECK_LT( maximumSeparationAberrationCorrection * radiansToArcseconds, 2.0e-2 );
    BOOST_CHECK_GT( maximumTransversePositionAngleAberrationCorrection * radiansToArcseconds, 1.0e-4 );
    BOOST_CHECK_LT( maximumTransversePositionAngleAberrationCorrection * radiansToArcseconds, 1.0e-3 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat

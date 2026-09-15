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
#include <sstream>
#include <string>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/matrixTextFileReader.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"

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
    }

    const double radiansToMilliarcseconds = 180.0 / mathematical_constants::PI * 3600.0 * 1000.0;
    const double separationResidualRms =
            std::sqrt( squaredSeparationResidualSum / static_cast< double >( observations.rows( ) ) ) * radiansToMilliarcseconds;
    const double transversePositionAngleResidualRms =
            std::sqrt( squaredTransversePositionAngleResidualSum / static_cast< double >( observations.rows( ) ) ) *
            radiansToMilliarcseconds;

    // PLU060 gives milliarcsecond-level pre-fit residuals comparable to the 1997 fitted-orbit
    // residuals (2.93 mas radial and 2.70 mas transverse) stored in the source data file.
    BOOST_CHECK_SMALL( separationResidualRms - 4.2227973382, 1.0e-6 );
    BOOST_CHECK_SMALL( transversePositionAngleResidualRms - 4.4095016757, 1.0e-6 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat

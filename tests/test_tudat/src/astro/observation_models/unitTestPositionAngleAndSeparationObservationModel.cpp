/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"
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
    std::shared_ptr< ObservationModelSettings > paSettings =
            std::make_shared< PositionAngleObservationModelSettings >( linkEnds, lightTimeCorrectionSettings );
    std::shared_ptr< ObservationModelSettings > sepSettings =
            std::make_shared< SeparationObservationModelSettings >( linkEnds, lightTimeCorrectionSettings );
    std::shared_ptr< ObservationModelSettings > pasSettings =
            std::make_shared< PositionAngleAndSeparationObservationModelSettings >( linkEnds, lightTimeCorrectionSettings );

    // Create observation models
    std::shared_ptr< ObservationModel< 1 > > paModel =
            ObservationModelCreator< 1, double, double >::createObservationModel( paSettings, bodies );
    std::shared_ptr< ObservationModel< 1 > > sepModel =
            ObservationModelCreator< 1, double, double >::createObservationModel( sepSettings, bodies );
    std::shared_ptr< ObservationModel< 2 > > pasModel =
            ObservationModelCreator< 2, double, double >::createObservationModel( pasSettings, bodies );

    // Test at several epochs
    double receiverObservationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;

    std::vector< double > linkEndTimesPA, linkEndTimesSep, linkEndTimesPAS;
    std::vector< Eigen::Vector6d > linkEndStatesPA, linkEndStatesSep, linkEndStatesPAS;

    Eigen::VectorXd positionAngleObservation =
            paModel->computeObservationsWithLinkEndData( receiverObservationTime, receiver, linkEndTimesPA, linkEndStatesPA );
    Eigen::VectorXd separationObservation =
            sepModel->computeObservationsWithLinkEndData( receiverObservationTime, receiver, linkEndTimesSep, linkEndStatesSep );
    Eigen::VectorXd positionAngleAndSeparationObservation =
            pasModel->computeObservationsWithLinkEndData( receiverObservationTime, receiver, linkEndTimesPAS, linkEndStatesPAS );

    // Verify combined model matches individual models
    BOOST_CHECK_CLOSE_FRACTION(
            positionAngleObservation( 0 ), positionAngleAndSeparationObservation( 0 ), std::numeric_limits< double >::epsilon( ) * 10.0 );
    BOOST_CHECK_CLOSE_FRACTION(
            separationObservation( 0 ), positionAngleAndSeparationObservation( 1 ), std::numeric_limits< double >::epsilon( ) * 10.0 );

    // Verify link end states match
    for( int i = 0; i < 3; i++ )
    {
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( linkEndStatesPA[ i ], linkEndStatesPAS[ i ], std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( linkEndStatesSep[ i ], linkEndStatesPAS[ i ], std::numeric_limits< double >::epsilon( ) );
    }

    // ===== Cross-consistency check with RA/DEC (relative angular position) model =====
    // Compute RA/DEC from the relative angular position model's link end states,
    // then convert to position angle and separation using spherical trig.
    std::shared_ptr< ObservationModelSettings > relAngPosSettings =
            std::make_shared< ObservationModelSettings >( relative_angular_position, linkEnds, lightTimeCorrectionSettings );
    std::shared_ptr< ObservationModel< 2 > > relAngPosModel =
            ObservationModelCreator< 2, double, double >::createObservationModel( relAngPosSettings, bodies );

    std::vector< double > linkEndTimesRelAngPos;
    std::vector< Eigen::Vector6d > linkEndStatesRelAngPos;
    relAngPosModel->computeObservationsWithLinkEndData( receiverObservationTime, receiver, linkEndTimesRelAngPos, linkEndStatesRelAngPos );

    // linkEndStatesRelAngPos: [0]=Mars (transmitter), [1]=Phobos (transmitter2), [2]=Earth (receiver)
    auto computeRaDec = []( const Eigen::Vector6d& transmitterState, const Eigen::Vector6d& receiverState ) {
        Eigen::Vector3d posDiff = ( transmitterState - receiverState ).segment( 0, 3 );
        double ra =
                2.0 * std::atan( posDiff( 1 ) / ( std::sqrt( posDiff( 0 ) * posDiff( 0 ) + posDiff( 1 ) * posDiff( 1 ) ) + posDiff( 0 ) ) );
        double dec = mathematical_constants::PI / 2.0 - std::acos( posDiff( 2 ) / posDiff.norm( ) );
        return std::make_pair( ra, dec );
    };

    auto [ raMars, decMars ] = computeRaDec( linkEndStatesRelAngPos[ 0 ], linkEndStatesRelAngPos[ 2 ] );
    auto [ raPhobos, decPhobos ] = computeRaDec( linkEndStatesRelAngPos[ 1 ], linkEndStatesRelAngPos[ 2 ] );

    // Convert RA/DEC to position angle and separation
    double deltaRa = raPhobos - raMars;
    double computedSeparation =
            std::acos( std::sin( decMars ) * std::sin( decPhobos ) + std::cos( decMars ) * std::cos( decPhobos ) * std::cos( deltaRa ) );
    double computedPositionAngle =
            std::atan2( std::cos( decPhobos ) * std::sin( deltaRa ),
                        std::sin( decPhobos ) * std::cos( decMars ) - std::cos( decPhobos ) * std::sin( decMars ) * std::cos( deltaRa ) );

    // Compare against dedicated models
    BOOST_CHECK_SMALL( computedPositionAngle - positionAngleObservation( 0 ), 1.0e-11 );
    BOOST_CHECK_SMALL( computedSeparation - separationObservation( 0 ), 1.0e-11 );
    BOOST_CHECK_SMALL( computedPositionAngle - positionAngleAndSeparationObservation( 0 ), 1.0e-11 );
    BOOST_CHECK_SMALL( computedSeparation - positionAngleAndSeparationObservation( 1 ), 1.0e-11 );

    // Test error: wrong reference link end
    BOOST_CHECK_THROW( paModel->computeObservationsWithLinkEndData( receiverObservationTime, transmitter, linkEndTimesPA, linkEndStatesPA ),
                       std::runtime_error );
    BOOST_CHECK_THROW(
            pasModel->computeObservationsWithLinkEndData( receiverObservationTime, transmitter, linkEndTimesPAS, linkEndStatesPAS ),
            std::runtime_error );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat

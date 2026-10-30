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

#include <limits>
#include <string>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationRate.h"
#include "tudat/simulation/estimation_setup/createObservationPartials.h"
#include "tudat/support/numericalObservationPartial.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/support/observationPartialTestFunctions.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/observation_models/stellarAberrationCorrection.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::observation_partials;
using namespace tudat::estimatable_parameters;

BOOST_AUTO_TEST_SUITE( test_position_angle_and_separation_partials )

//! Check the Jacobians of the stellar-aberration transformation against central differences.
BOOST_AUTO_TEST_CASE( testStellarAberrationPartials )
{
    // A large velocity (about 0.01c) makes the aberration factors clearly distinct from the identity.
    const Eigen::Vector3d astrometricDirection = Eigen::Vector3d( 3.0, -1.0, 2.0 ) * 1.0E11;
    const Eigen::Vector3d observerVelocity( 1.2E6, -2.1E6, 0.8E6 );

    Eigen::Matrix3d partialWrtAstrometricDirection;
    Eigen::Matrix3d partialWrtObserverVelocity;
    calculateStellarAberrationPartials(
            astrometricDirection, observerVelocity, partialWrtAstrometricDirection, partialWrtObserverVelocity );

    Eigen::Matrix3d numericalPartialWrtAstrometricDirection;
    Eigen::Matrix3d numericalPartialWrtObserverVelocity;
    const double directionPerturbation = 1.0E5;
    const double velocityPerturbation = 1.0E2;
    for( int i = 0; i < 3; i++ )
    {
        Eigen::Vector3d directionChange = Eigen::Vector3d::Zero( );
        directionChange( i ) = directionPerturbation;
        numericalPartialWrtAstrometricDirection.col( i ) =
                ( calculateApparentDirectionWithStellarAberration( astrometricDirection + directionChange, observerVelocity ) -
                  calculateApparentDirectionWithStellarAberration( astrometricDirection - directionChange, observerVelocity ) ) /
                ( 2.0 * directionPerturbation );

        Eigen::Vector3d velocityChange = Eigen::Vector3d::Zero( );
        velocityChange( i ) = velocityPerturbation;
        numericalPartialWrtObserverVelocity.col( i ) =
                ( calculateApparentDirectionWithStellarAberration( astrometricDirection, observerVelocity + velocityChange ) -
                  calculateApparentDirectionWithStellarAberration( astrometricDirection, observerVelocity - velocityChange ) ) /
                ( 2.0 * velocityPerturbation );
    }
    checkMatrixCloseRelativeToNorm( partialWrtAstrometricDirection, numericalPartialWrtAstrometricDirection, 1.0E-8 );
    checkMatrixCloseRelativeToNorm( partialWrtObserverVelocity, numericalPartialWrtObserverVelocity, 1.0E-8 );

    // The aberration factors must differ from those of a static observer by much more than the tolerance above.
    Eigen::Matrix3d staticObserverPartialWrtAstrometricDirection;
    Eigen::Matrix3d staticObserverPartialWrtObserverVelocity;
    calculateStellarAberrationPartials( astrometricDirection,
                                        Eigen::Vector3d::Zero( ),
                                        staticObserverPartialWrtAstrometricDirection,
                                        staticObserverPartialWrtObserverVelocity );
    BOOST_CHECK_GT( ( partialWrtAstrometricDirection - staticObserverPartialWrtAstrometricDirection ).norm( ) /
                            staticObserverPartialWrtAstrometricDirection.norm( ),
                    1.0E-3 );
}

//! Compare analytical and numerical partials of aberrated P/S observables, including the receiver-velocity contribution.
BOOST_AUTO_TEST_CASE( testAberratedPositionAngleAndSeparationPartials )
{
    LinkDefinition linkEnds;
    linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "Graz" );
    linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "MSL" );
    linkEnds[ transmitter2 ] = std::make_pair< std::string, std::string >( "Moon", "" );

    const std::vector< std::shared_ptr< ObservationAncillarySimulationSettings > > ancillarySettingsList = {
        getPositionAngleAncillarySettings( j2000_position_angle_reference_frame, TUDAT_NAN, true ),
        getPositionAngleAncillarySettings( true_of_date_iau_2006_2000a_position_angle_reference_frame, TUDAT_NAN, true )
    };
    for( const auto& ancillarySettings : ancillarySettingsList )
    {
        testAberratedObservationPartials< 1 >( std::make_shared< ObservationModelSettings >( position_angle, linkEnds ),
                                               ancillarySettings );
        testAberratedObservationPartials< 1 >( std::make_shared< ObservationModelSettings >( separation_distance, linkEnds ),
                                               ancillarySettings );
        testAberratedObservationPartials< 2 >( std::make_shared< ObservationModelSettings >( position_angle_and_separation, linkEnds ),
                                               ancillarySettings );
    }
}

BOOST_AUTO_TEST_CASE( testPositionAngleAndSeparationPartials )
{
    // Define and create ground stations.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.resize( 3 );
    groundStations[ 0 ] = std::make_pair( "Earth", "Graz" );
    groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );
    groundStations[ 2 ] = std::make_pair( "Moon", "" );

    std::vector< std::string > perturbingBodies;
    perturbingBodies.push_back( "Earth" );

    // Exercise the same numerical-partial comparisons for each celestial pole.
    // The custom pole is intentionally tilted by more than 50 degrees from J2000 north.
    const std::vector< std::pair< std::string, std::shared_ptr< ObservationAncillarySimulationSettings > > > referencePoleCases = {
        { "default J2000", nullptr },
        { "B1950", getPositionAngleAncillarySettings( b1950_position_angle_reference_frame ) },
        { "true of date", getPositionAngleAncillarySettings( true_of_date_iau_1976_1980_position_angle_reference_frame ) },
        { "custom", getCustomPositionAngleAncillarySettings( Eigen::Vector3d( 0.4, -0.7, 0.6 ).normalized( ) ) }
    };
    for( const auto& referencePoleCase : referencePoleCases )
    {
        BOOST_TEST_MESSAGE( "Testing P/S partials with " << referencePoleCase.first << " reference pole" );

        // Test partials with constant ephemerides (allows test of position partials)
        {
            // Create environment
            SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true );

            // Set link ends for observation model
            LinkDefinition linkEnds;
            linkEnds[ receiver ] = groundStations[ 0 ];
            linkEnds[ transmitter ] = groundStations[ 1 ];
            linkEnds[ transmitter2 ] = groundStations[ 2 ];

            // Test position angle partials
            {
                std::shared_ptr< ObservationModel< 1 > > positionAngleModel =
                        observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                                std::make_shared< observation_models::ObservationModelSettings >(
                                        observation_models::position_angle,
                                        linkEnds,
                                        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >{
                                                std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                                        perturbingBodies ) } ),
                                bodies );

                std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                        createEstimatableParameters( bodies, 1.1E7 );

                testObservationPartials< 1 >( positionAngleModel,
                                              bodies,
                                              fullEstimatableParameterSet,
                                              linkEnds,
                                              position_angle,
                                              1.0E-4,
                                              true,
                                              true,
                                              1.0,
                                              Eigen::VectorXd::Constant( 4, 1.0 ),
                                              referencePoleCase.second );
            }

            // Test separation distance partials
            {
                std::shared_ptr< ObservationModel< 1 > > separationModel =
                        observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                                std::make_shared< observation_models::ObservationModelSettings >(
                                        observation_models::separation_distance,
                                        linkEnds,
                                        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >{
                                                std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                                        perturbingBodies ) } ),
                                bodies );

                std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                        createEstimatableParameters( bodies, 1.1E7 );

                testObservationPartials< 1 >( separationModel,
                                              bodies,
                                              fullEstimatableParameterSet,
                                              linkEnds,
                                              separation_distance,
                                              1.0E-4,
                                              true,
                                              true,
                                              1.0,
                                              Eigen::VectorXd::Constant( 4, 1.0 ),
                                              referencePoleCase.second );
            }

            // Test combined position angle and separation distance partials
            {
                std::shared_ptr< ObservationModel< 2 > > positionAngleAndSeparationModel =
                        observation_models::ObservationModelCreator< 2, double, double >::createObservationModel(
                                std::make_shared< observation_models::ObservationModelSettings >(
                                        observation_models::position_angle_and_separation,
                                        linkEnds,
                                        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >{
                                                std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                                        perturbingBodies ) } ),
                                bodies );

                std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                        createEstimatableParameters( bodies, 1.1E7 );

                testObservationPartials< 2 >( positionAngleAndSeparationModel,
                                              bodies,
                                              fullEstimatableParameterSet,
                                              linkEnds,
                                              position_angle_and_separation,
                                              1.0E-4,
                                              true,
                                              true,
                                              1.0,
                                              Eigen::VectorXd::Constant( 4, 1.0 ),
                                              referencePoleCase.second );
            }
        }

        // Test partials with real ephemerides (without test of position partials)
        {
            // Create environment
            SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, false );

            // Set link ends for observation model
            LinkDefinition linkEnds;
            linkEnds[ receiver ] = groundStations[ 0 ];
            linkEnds[ transmitter ] = groundStations[ 1 ];
            linkEnds[ transmitter2 ] = groundStations[ 2 ];

            // Test position angle partials
            {
                std::shared_ptr< ObservationModel< 1 > > positionAngleModel =
                        observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                                std::make_shared< observation_models::ObservationModelSettings >(
                                        observation_models::position_angle,
                                        linkEnds,
                                        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >{
                                                std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                                        perturbingBodies ) } ),
                                bodies );

                std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                        createEstimatableParameters( bodies, 1.1E7 );

                testObservationPartials< 1 >( positionAngleModel,
                                              bodies,
                                              fullEstimatableParameterSet,
                                              linkEnds,
                                              position_angle,
                                              1.0E-4,
                                              false,
                                              true,
                                              1.0,
                                              Eigen::VectorXd::Constant( 4, 1.0 ),
                                              referencePoleCase.second );
            }

            // Test separation distance partials
            {
                std::shared_ptr< ObservationModel< 1 > > separationModel =
                        observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                                std::make_shared< observation_models::ObservationModelSettings >(
                                        observation_models::separation_distance,
                                        linkEnds,
                                        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >{
                                                std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                                        perturbingBodies ) } ),
                                bodies );

                std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                        createEstimatableParameters( bodies, 1.1E7 );

                testObservationPartials< 1 >( separationModel,
                                              bodies,
                                              fullEstimatableParameterSet,
                                              linkEnds,
                                              separation_distance,
                                              1.0E-4,
                                              false,
                                              true,
                                              1.0,
                                              Eigen::VectorXd::Constant( 4, 1.0 ),
                                              referencePoleCase.second );
            }

            // Test combined position angle and separation distance partials
            {
                std::shared_ptr< ObservationModel< 2 > > positionAngleAndSeparationModel =
                        observation_models::ObservationModelCreator< 2, double, double >::createObservationModel(
                                std::make_shared< observation_models::ObservationModelSettings >(
                                        observation_models::position_angle_and_separation,
                                        linkEnds,
                                        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >{
                                                std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                                        perturbingBodies ) } ),
                                bodies );

                std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                        createEstimatableParameters( bodies, 1.1E7 );

                testObservationPartials< 2 >( positionAngleAndSeparationModel,
                                              bodies,
                                              fullEstimatableParameterSet,
                                              linkEnds,
                                              position_angle_and_separation,
                                              1.0E-4,
                                              false,
                                              true,
                                              1.0,
                                              Eigen::VectorXd::Constant( 4, 1.0 ),
                                              referencePoleCase.second );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat

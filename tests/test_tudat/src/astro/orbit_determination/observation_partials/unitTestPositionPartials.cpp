/*    Copyright (c) 2010-2019, Delft University of Technology
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
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationRate.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationalOrientation.h"
#include "tudat/simulation/estimation_setup/createCartesianStatePartials.h"
#include "tudat/support/numericalObservationPartial.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/environment_setup/createBodyDeformationModel.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/displacementLoveNumbers.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::basic_astrodynamics;
using namespace tudat::ground_stations;
using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::observation_partials;

BOOST_AUTO_TEST_SUITE( test_position_partials )

//! Test partial derivatives of positions w.r.t. parameters. Partials of most observables are computed in terms of these
//! partials
BOOST_AUTO_TEST_CASE( testCartesianStatePartials )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Specify initial time
    double initialEphemerisTime = 1.0E7;

    // Create bodies.
    SystemOfBodies bodies;
    bodies.createEmptyBody( "Earth" );
    bodies.createEmptyBody( "Moon" );
    bodies.createEmptyBody( "Sun" );
    bodies.createEmptyBody( "Mars" );

    // Define properties of bodies
    bodies.at( "Earth" )->setShapeModel( std::make_shared< SphericalBodyShapeModel >( spice_interface::getAverageRadius( "Earth" ) ) );
    bodies.at( "Mars" )->setShapeModel( std::make_shared< SphericalBodyShapeModel >( spice_interface::getAverageRadius( "Mars" ) ) );

    bodies.at( "Earth" )->setEphemeris( std::make_shared< ConstantEphemeris >(
            getBodyCartesianStateAtEpoch( "Earth", "SSB", "ECLIPJ2000", "NONE", 0.0 ), "SSB", "ECLIPJ2000" ) );
    bodies.at( "Moon" )->setEphemeris( std::make_shared< ConstantEphemeris >(
            getBodyCartesianStateAtEpoch( "Moon", "SSB", "ECLIPJ2000", "NONE", 0.0 ), "SSB", "ECLIPJ2000" ) );
    bodies.at( "Sun" )->setEphemeris( std::make_shared< ConstantEphemeris >(
            getBodyCartesianStateAtEpoch( "Sun", "SSB", "ECLIPJ2000", "NONE", 0.0 ), "SSB", "ECLIPJ2000" ) );
    bodies.at( "Mars" )->setEphemeris( std::make_shared< ConstantEphemeris >(
            getBodyCartesianStateAtEpoch( "Mars", "SSB", "ECLIPJ2000", "NONE", 0.0 ), "SSB", "ECLIPJ2000" ) );

    bodies.at( "Sun" )->setGravityFieldModel( std::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Sun" ) ) );
    bodies.at( "Moon" )->setGravityFieldModel( std::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Moon" ) ) );
    bodies.at( "Earth" )->setGravityFieldModel( std::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Earth" ) ) );

    bodies.at( "Earth" )->setRotationalEphemeris( createRotationModel(
            std::make_shared< SimpleRotationModelSettings >(
                    "ECLIPJ2000",
                    "IAU_Earth",
                    spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
                    initialEphemerisTime,
                    2.0 * mathematical_constants::PI / physical_constants::JULIAN_DAY ),
            "Earth" ) );
    bodies.at( "Mars" )->setRotationalEphemeris( createRotationModel(
            std::make_shared< SimpleRotationModelSettings >(
                    "ECLIPJ2000",
                    "IAU_Mars",
                    spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "IAU_Mars", initialEphemerisTime ),
                    initialEphemerisTime,
                    2.0 * mathematical_constants::PI / physical_constants::JULIAN_DAY ),
            "Mars" ) );

    // Create ground stations
    std::map< std::pair< std::string, std::string >, Eigen::Vector3d > groundStationsToCreate;
    groundStationsToCreate[ std::make_pair( "Earth", "Graz" ) ] = ( Eigen::Vector3d( ) << 1.7E6, -6.2E6, 1.3E5 ).finished( );
    groundStationsToCreate[ std::make_pair( "Mars", "MSL" ) ] = ( Eigen::Vector3d( ) << -2.5E5, 3.2E6, -2.65E4 ).finished( );
    createGroundStations( bodies, groundStationsToCreate );

    // Define list of ground station names.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.push_back( std::make_pair( "Earth", "Graz" ) );
    groundStations.push_back( std::make_pair( "Mars", "MSL" ) );

    // Create link ends set.
    LinkDefinition linkEnds;
    linkEnds[ observed_body ] = groundStations[ 0 ];

    LinkDefinition linkEnds2;
    linkEnds2[ observed_body ] = groundStations[ 1 ];

    std::shared_ptr< GroundStation > receivingGroundStation = bodies.at( "Earth" )->getGroundStation( "Graz" );

    // Create parameter objects.
    std::shared_ptr< RotationRate > earthRotationRate = std::make_shared< RotationRate >(
            std::dynamic_pointer_cast< SimpleRotationalEphemeris >( bodies.at( "Earth" )->getRotationalEphemeris( ) ), "Earth" );
    std::shared_ptr< ConstantRotationalOrientation > earthPolePosition = std::make_shared< ConstantRotationalOrientation >(
            std::dynamic_pointer_cast< SimpleRotationalEphemeris >( bodies.at( "Earth" )->getRotationalEphemeris( ) ), "Earth" );

    // Create explicit position partial objects.
    std::shared_ptr< CartesianStatePartial > partialObjectWrtReceiverPosition =
            createCartesianStatePartialsWrtBodyState( linkEnds.linkEnds_, bodies, "Earth" ).begin( )->second;

    // Create explicit parameter partial objects.
    std::shared_ptr< CartesianStatePartial > partialObjectWrtReceiverRotationRate =
            createCartesianStatePartialsWrtParameter( linkEnds.linkEnds_, bodies, earthRotationRate ).begin( )->second;
    std::shared_ptr< CartesianStatePartial > partialObjectWrtReceiverPolePosition =
            createCartesianStatePartialsWrtParameter( linkEnds.linkEnds_, bodies, earthPolePosition ).begin( )->second;

    // Calculate transmission/reception times and states
    Eigen::Vector6d currentState;
    double receptionTime = 1.1E7;
    currentState = bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris( receptionTime );

    double currentTime = receptionTime;

    // Compute partials
    Eigen::MatrixXd partialWrtReceiverPosition = partialObjectWrtReceiverPosition->calculatePartialOfPosition( currentState, currentTime );

    Eigen::MatrixXd partialWrtReceiverRotationRate =
            partialObjectWrtReceiverRotationRate->calculatePartialOfPosition( currentState, currentTime );
    Eigen::MatrixXd partialOfVelocityWrtReceiverRotationRate =
            partialObjectWrtReceiverRotationRate->calculatePartialOfVelocity( currentState, currentTime );

    Eigen::MatrixXd partialWrtReceiverPolePosition =
            partialObjectWrtReceiverPolePosition->calculatePartialOfPosition( currentState, currentTime );
    Eigen::MatrixXd partialOfVelocityWrtReceiverPolePosition =
            partialObjectWrtReceiverPolePosition->calculatePartialOfVelocity( currentState, currentTime );

    // Define observation function
    std::function< Eigen::VectorXd( const double ) > observationFunctionAtReception =
            std::bind( &Ephemeris::getCartesianState,
                       createReferencePointCompositeEphemeris< double, double >(
                               bodies.at( "Earth" )->getEphemeris( ),
                               bodies.at( "Earth" )->getRotationalEphemeris( ),
                               std::bind( &GroundStation::getStateInPlanetFixedFrame< double, double >,
                                          bodies.at( "Earth" )->getGroundStation( "Graz" ),
                                          std::placeholders::_1,
                                          "Earth" ) ),
                       std::placeholders::_1 );

    // Calculate numerical partials w.r.t. Earth state.
    Eigen::Vector3d bodyPositionVariation;
    bodyPositionVariation << 10.0, 10.0, 10.0;
    std::shared_ptr< ConstantEphemeris > earthEphemeris =
            std::dynamic_pointer_cast< ConstantEphemeris >( bodies.at( "Earth" )->getEphemeris( ) );
    Eigen::Vector6d earthUnperturbedState = earthEphemeris->getCartesianState( 0.0 );
    Eigen::Vector6d perturbedEarthState;
    Eigen::Matrix< double, 3, 3 > numericalPartialWrtReceiverPosition = Eigen::Matrix< double, 3, 3 >::Zero( );
    for( int i = 0; i < 3; i++ )
    {
        perturbedEarthState = earthUnperturbedState;
        perturbedEarthState( i ) += bodyPositionVariation( i );
        earthEphemeris->updateConstantState( perturbedEarthState );
        Eigen::Vector3d upPerturbedPosition = observationFunctionAtReception( currentTime ).segment( 0, 3 );

        perturbedEarthState = earthUnperturbedState;
        perturbedEarthState( i ) -= bodyPositionVariation( i );
        earthEphemeris->updateConstantState( perturbedEarthState );
        Eigen::Vector3d downPerturbedPosition = observationFunctionAtReception( currentTime ).segment( 0, 3 );

        numericalPartialWrtReceiverPosition.block( 0, i, 3, 1 ) =
                ( upPerturbedPosition - downPerturbedPosition ) / ( 2.0 * bodyPositionVariation( i ) );
    }
    earthEphemeris->updateConstantState( earthUnperturbedState );

    // Test partial w.r.t. position
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtReceiverPosition.block( 0, 0, 3, 3 ), numericalPartialWrtReceiverPosition, 1.0E-12 );

    // Compute numerical partial w.r.t. rotation rate.
    Eigen::Vector6d numericalPartialWrtReceiverRotationRate =
            calculateNumericalObservationParameterPartial( earthRotationRate, 1.0E-10, observationFunctionAtReception, receptionTime );

    // Test partial w.r.t. rotation rate
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtReceiverRotationRate, numericalPartialWrtReceiverRotationRate.segment( 0, 3 ), 1.0E-5 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            partialOfVelocityWrtReceiverRotationRate, numericalPartialWrtReceiverRotationRate.segment( 3, 3 ), 1.0E-5 );

    // Compute numerical partial w.r.t. pole position
    Eigen::VectorXd polePositionPerturbation = ( Eigen::Vector2d( ) << 1.0E-5, 1.0E-5 ).finished( );
    Eigen::MatrixXd numericalPartialWrtReceiverPolePosition = calculateNumericalObservationParameterPartial(
            earthPolePosition, polePositionPerturbation, observationFunctionAtReception, receptionTime );

    // Test partials w.r.t. pole position (different tolernaces due to different magnitudes of partials).
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( partialWrtReceiverPolePosition.block( 0, 0, 1, 2 ) ),
                                       ( numericalPartialWrtReceiverPolePosition.block( 0, 0, 1, 2 ) ),
                                       1.0E-4 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( partialWrtReceiverPolePosition.block( 1, 0, 2, 2 ) ),
                                       ( numericalPartialWrtReceiverPolePosition.block( 1, 0, 2, 2 ) ),
                                       1.0E-6 );

    for( int i = 0; i < 3; i++ )
    {
        for( int j = 0; j < 2; j++ )
        {
            BOOST_CHECK_SMALL(
                    std::fabs( partialOfVelocityWrtReceiverPolePosition( i, j ) - numericalPartialWrtReceiverPolePosition( i + 3, j ) ),
                    1.0E-6 );
        }
    }
}

//! Test partial derivatives of the position and velocity of a ground station w.r.t. the displacement Love and Shida numbers of the
//! basic tidal deformation of its body
BOOST_AUTO_TEST_CASE( testDisplacementLoveNumberPartials )
{
    spice_interface::loadStandardSpiceKernels( );

    double initialEphemerisTime = 1.0E7;

    // Create bodies, with a basic tidal deformation of the Earth due to the Moon and the Sun. The (constant) states are relative to the
    // Earth, to limit rounding errors in the numerical partials.
    SystemOfBodies bodies;
    bodies.createEmptyBody( "Earth" );
    bodies.createEmptyBody( "Moon" );
    bodies.createEmptyBody( "Sun" );
    for( const std::string bodyName : { "Earth", "Moon", "Sun" } )
    {
        bodies.at( bodyName )
                ->setEphemeris( std::make_shared< ConstantEphemeris >(
                        getBodyCartesianStateAtEpoch( bodyName, "Earth", "ECLIPJ2000", "NONE", 0.0 ), "SSB", "ECLIPJ2000" ) );
        bodies.at( bodyName )->setGravityFieldModel( std::make_shared< GravityFieldModel >( getBodyGravitationalParameter( bodyName ) ) );
    }
    bodies.at( "Earth" )->setRotationalEphemeris( createRotationModel(
            std::make_shared< SimpleRotationModelSettings >(
                    "ECLIPJ2000",
                    "IAU_Earth",
                    spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
                    initialEphemerisTime,
                    2.0 * mathematical_constants::PI / physical_constants::JULIAN_DAY ),
            "Earth" ) );
    std::map< int, std::pair< double, double > > displacementLoveNumbers;
    displacementLoveNumbers[ 2 ] = std::make_pair( 0.6, 0.08 );
    displacementLoveNumbers[ 3 ] = std::make_pair( 0.29, 0.015 );
    bodies.at( "Earth" )->addBodyDeformationModel( createBodyDeformationModel(
            basicTidalBodyShapeDeformation( { "Moon", "Sun" }, displacementLoveNumbers, 6378137.0 ), "Earth", bodies ) );

    // Create ground station, of which the position includes the tidal displacement
    createGroundStation( bodies.at( "Earth" ),
                         "Graz",
                         ( Eigen::Vector3d( ) << 1.7E6, -6.2E6, 1.3E5 ).finished( ),
                         coordinate_conversions::cartesian_position,
                         { bodyDeformationStationMotionSettings( ) } );
    LinkDefinition linkEnds;
    linkEnds[ observed_body ] = std::make_pair< std::string, std::string >( "Earth", "Graz" );

    std::function< Eigen::VectorXd( const double ) > observationFunctionAtReception =
            std::bind( &Ephemeris::getCartesianState,
                       createReferencePointCompositeEphemeris< double, double >(
                               bodies.at( "Earth" )->getEphemeris( ),
                               bodies.at( "Earth" )->getRotationalEphemeris( ),
                               std::bind( &GroundStation::getStateInPlanetFixedFrame< double, double >,
                                          bodies.at( "Earth" )->getGroundStation( "Graz" ),
                                          std::placeholders::_1,
                                          "Earth" ) ),
                       std::placeholders::_1 );

    std::vector< std::shared_ptr< BasicTidalBodyDeformation > > deformationModels = {
        std::dynamic_pointer_cast< BasicTidalBodyDeformation >( bodies.at( "Earth" )->getBodyDeformationModels( ).at( 0 ) )
    };
    double receptionTime = 1.1E7;
    Eigen::Vector6d currentState = bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris( receptionTime );

    // Compare analytical and numerical partials w.r.t. the Love numbers h_2, h_3 and the Shida numbers l_2, l_3
    for( const bool useShidaNumbers : { false, true } )
    {
        std::shared_ptr< DisplacementLoveNumbers > loveNumberParameter =
                std::make_shared< DisplacementLoveNumbers >( deformationModels, std::vector< int >( { 2, 3 } ), useShidaNumbers, "Earth" );
        BOOST_CHECK_EQUAL( loveNumberParameter->getParameterSize( ), 2 );

        std::shared_ptr< CartesianStatePartial > partialObject =
                createCartesianStatePartialsWrtParameter( linkEnds.linkEnds_, bodies, loveNumberParameter ).begin( )->second;
        Eigen::MatrixXd partialOfPosition = partialObject->calculatePartialOfPosition( currentState, receptionTime );
        Eigen::MatrixXd partialOfVelocity = partialObject->calculatePartialOfVelocity( currentState, receptionTime );

        // The displacement is linear in the Love and Shida numbers, so that a large perturbation can be used
        Eigen::VectorXd perturbation = Eigen::VectorXd::Constant( 2, 1.0 );
        Eigen::MatrixXd numericalPartial = calculateNumericalObservationParameterPartial(
                loveNumberParameter, perturbation, observationFunctionAtReception, receptionTime );

        // Compare per degree, relative to the size of the partial (millimetres to decimetres per unit Love or Shida number)
        for( int i = 0; i < 2; i++ )
        {
            BOOST_CHECK( partialOfPosition.col( i ).norm( ) > 1.0E-3 );
            BOOST_CHECK_SMALL(
                    ( partialOfPosition.col( i ) - numericalPartial.block( 0, i, 3, 1 ) ).norm( ) / partialOfPosition.col( i ).norm( ),
                    1.0E-6 );
            BOOST_CHECK_SMALL(
                    ( partialOfVelocity.col( i ) - numericalPartial.block( 3, i, 3, 1 ) ).norm( ) / partialOfVelocity.col( i ).norm( ),
                    1.0E-6 );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat

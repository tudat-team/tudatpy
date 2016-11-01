
/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *
 *    References
 *
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"
#include "Tudat/Astrodynamics/Ephemerides/constantEphemeris.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/constantRotationRate.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/constantRotationalOrientation.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/createPositionPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/UnitTests/numericalObservationPartial.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/EstimationSetup/createEstimatableParameters.h"


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

BOOST_AUTO_TEST_SUITE( test_position_partials)

BOOST_AUTO_TEST_CASE( testPositionPartials )
{

    //Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    loadSpiceKernelInTudat( kernelsPath + "de421.bsp");
    loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");

    // Specify initial time
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = 1.2E7;
    double buffer = 1000.0;

    // Create bodies.
    NamedBodyMap bodyMap;
    bodyMap[ "Earth" ] = boost::make_shared< Body >( );
    bodyMap[ "Moon" ] = boost::make_shared< Body >( );
    bodyMap[ "Sun" ] = boost::make_shared< Body >( );
    bodyMap[ "Mars" ] = boost::make_shared< Body >( );

    ( bodyMap[ "Earth" ] )->setShapeModel( boost::make_shared< SphericalBodyShapeModel >( spice_interface::getAverageRadius( "Earth" ) ) );
    ( bodyMap[ "Mars" ] )->setShapeModel( boost::make_shared< SphericalBodyShapeModel >( spice_interface::getAverageRadius( "Mars" ) ) );

    bodyMap[ "Earth" ]->setEphemeris(
                boost::make_shared< ConstantEphemeris >(
                    getBodyCartesianStateAtEpoch( "Earth", "SSB", "ECLIPJ2000", "NONE", 0.0 ),
                    "SSB", "ECLIPJ2000" ) );
    bodyMap[ "Moon" ]->setEphemeris(
                boost::make_shared< ConstantEphemeris >(
                    getBodyCartesianStateAtEpoch( "Moon", "SSB", "ECLIPJ2000", "NONE", 0.0 ),
                    "SSB", "ECLIPJ2000" ) );
    bodyMap[ "Sun" ]->setEphemeris(
                boost::make_shared< ConstantEphemeris >(
                    getBodyCartesianStateAtEpoch( "Sun", "SSB", "ECLIPJ2000", "NONE", 0.0 ),
                    "SSB", "ECLIPJ2000" ) );
    bodyMap[ "Mars" ]->setEphemeris(
                boost::make_shared< ConstantEphemeris >(
                    getBodyCartesianStateAtEpoch( "Mars", "SSB", "ECLIPJ2000", "NONE", 0.0 ),
                    "SSB", "ECLIPJ2000" ) );

    ( bodyMap[ "Sun" ] )->setGravityFieldModel(
                boost::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Sun" ) ) );
    ( bodyMap[ "Moon" ] )->setGravityFieldModel(
                boost::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Moon" ) ) );
    ( bodyMap[ "Earth" ] )->setGravityFieldModel(
                boost::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Earth" ) ) );



    ( bodyMap[ "Earth" ] )->setRotationalEphemeris(
                createRotationModel(
                    boost::make_shared< SimpleRotationModelSettings >(
                        "ECLIPJ2000", "IAU_Earth",
                        spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
                        initialEphemerisTime, 2.0 * mathematical_constants::PI / physical_constants::JULIAN_DAY ), "Earth" ) );
    ( bodyMap[ "Mars" ] )->setRotationalEphemeris(
                createRotationModel(
                    boost::make_shared< SimpleRotationModelSettings >(
                        "ECLIPJ2000", "IAU_Mars",
                        spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "IAU_Mars", initialEphemerisTime ),
                        initialEphemerisTime, 2.0 * mathematical_constants::PI / physical_constants::JULIAN_DAY ), "Mars" ) );


    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


    std::map< std::pair< std::string, std::string >, Eigen::Vector3d > groundStationsToCreate;
    groundStationsToCreate[ std::make_pair( "Earth", "Graz" ) ] =
            ( Eigen::Vector3d( ) << 1.7E6, -6.2E6, 1.3E5 ).finished( );
    groundStationsToCreate[ std::make_pair( "Mars", "MSL" ) ] =
            ( Eigen::Vector3d( ) <<-2.5E5, 3.2E6, -2.65E4 ).finished( );


    createGroundStations( bodyMap, groundStationsToCreate );



    // Define and create ground stations.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.push_back( std::make_pair( "Earth", "Graz" ) );
    groundStations.push_back( std::make_pair( "Mars", "MSL" ) );


    LinkEnds linkEnds;
    linkEnds[ observed_body ] = groundStations[ 0 ];

    LinkEnds linkEnds2;
    linkEnds2[ observed_body ] = groundStations[ 1 ];

    boost::shared_ptr< GroundStation > receivingGroundStation =
            ( bodyMap[ "Earth" ] )->getGroundStation( "Graz" );

    //Create parameter objects.
    boost::shared_ptr< RotationRate > earthRotationRate = boost::make_shared< RotationRate >(
                boost::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                    bodyMap[ "Earth" ]->getRotationalEphemeris( ) ), "Earth");
    boost::shared_ptr< ConstantRotationalOrientation > earthPolePosition = boost::make_shared< ConstantRotationalOrientation >(
                boost::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                    bodyMap[ "Earth" ]->getRotationalEphemeris( ) ), "Earth" );




    // Create explicit range partial objects.
    boost::shared_ptr< PositionPartial > partialObjectWrtReceiverPosition =
            createPositionPartialsWrtBodyPosition( linkEnds, bodyMap, "Earth" ).begin( )->second;

    boost::shared_ptr< PositionPartial > partialObjectWrtReceiverRotationRate =
            createPositionPartialsWrtParameter(
                linkEnds, bodyMap, earthRotationRate ).begin( )->second;
    boost::shared_ptr< PositionPartial > partialObjectWrtReceiverPolePosition =
            createPositionPartialsWrtParameter(
                linkEnds, bodyMap, earthPolePosition ).begin( )->second;

    // Calculate transmission/reception times and states
    basic_mathematics::Vector6d currentState;
    double receptionTime = 1.1E7;
    currentState = bodyMap.at( "Earth" )->getStateInBaseFrameFromEphemeris( receptionTime );

    double currentTime = receptionTime;

    Eigen::MatrixXd partialWrtReceiverPosition =
            partialObjectWrtReceiverPosition->calculatePartial( currentState, currentTime );

    Eigen::MatrixXd partialWrtReceiverRotationRate =
            partialObjectWrtReceiverRotationRate->calculatePartial( currentState, currentTime );

    Eigen::MatrixXd partialWrtReceiverPolePosition =
            partialObjectWrtReceiverPolePosition->calculatePartial( currentState, currentTime );


    boost::function< Eigen::VectorXd( const double ) > observationFunctionAtReception =
            boost::bind( &Ephemeris::getCartesianState, createReferencePointEphemeris< double, double >(
                             bodyMap.at( "Earth" )->getEphemeris( ), bodyMap.at( "Earth" )->getRotationalEphemeris( ),
                             boost::bind( &GroundStation::getStateInPlanetFixedFrame< double, double >,
                                          bodyMap[ "Earth" ]->getGroundStation( "Graz" ), _1 ) ), _1 );


    double rotationRateVariation = 1.0E-10;

    Eigen::Vector3d bodyPositionVariation;
    bodyPositionVariation << 10.0, 10.0, 10.0;

    // Calculate numerical partials w.r.t. Earth state.
    boost::shared_ptr< ConstantEphemeris > earthEphemeris = boost::dynamic_pointer_cast< ConstantEphemeris >(
                bodyMap[ "Earth" ]->getEphemeris( ) );
    basic_mathematics::Vector6d earthUnperturbedState = earthEphemeris->getCartesianState( 0.0 );
    basic_mathematics::Vector6d perturbedEarthState;
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

        numericalPartialWrtReceiverPosition.block( 0, i, 3, 1 ) = ( upPerturbedPosition - downPerturbedPosition ) /
                ( 2.0 * bodyPositionVariation( i ) );
    }
    earthEphemeris->updateConstantState( earthUnperturbedState );

    Eigen::Vector3d numericalPartialWrtReceiverRotationRate = calculateNumericalObservationParameterPartial(
                earthRotationRate, 1.0E-10, observationFunctionAtReception,
                receptionTime );


    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtReceiverRotationRate, numericalPartialWrtReceiverRotationRate, 1.0E-5 );



    Eigen::VectorXd polePositionPerturbation = ( Eigen::Vector2d( )<<1.0E-5, 1.0E-5 ).finished( );
    Eigen::MatrixXd numericalPartialWrtReceiverPolePosition = calculateNumericalObservationParameterPartial(
                earthPolePosition, polePositionPerturbation, observationFunctionAtReception, receptionTime );


    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( partialWrtReceiverPolePosition.block( 0, 0, 1, 2 ) ),
                                       ( numericalPartialWrtReceiverPolePosition.block( 0, 0, 1, 2 ) ), 1.0E-4 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( partialWrtReceiverPolePosition.block( 1, 0, 2, 2 ) ),
                                       ( numericalPartialWrtReceiverPolePosition.block( 1, 0, 2, 2 ) ), 1.0E-6 );
}



BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat





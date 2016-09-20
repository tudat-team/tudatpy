
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

#include "Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"
#include "Astrodynamics/Ephemerides/constantEphemeris.h"
#include "Astrodynamics/Ephemerides/compositeRotationalEphemeris.h"
#include "Astrodynamics/Ephemerides/createLinkEndEphemeris.h"
#include "Astrodynamics/GroundStations/geodeticGroundStationState.h"
#include "Astrodynamics/OrbitDetermination/EstimatableParameters/displacementLoveNumbers.h"
#include "Astrodynamics/OrbitDetermination/EstimatableParameters/constantRotationRate.h"
#include "Astrodynamics/OrbitDetermination/EstimatableParameters/groundStationPosition.h"
#include "Astrodynamics/OrbitDetermination/EstimatableParameters/constantRotationalOrientation.h"
#include "Astrodynamics/OrbitDetermination/EstimatableParameters/moonLibrationParameters.h"
#include "Astrodynamics/OrbitDetermination/EstimatableParameters/rotationModelPeriodicVariationAmplitudes.h"
#include "Astrodynamics/OrbitDetermination/EstimatableParameters/rotationModelPolynomialVariations.h"
#include "Astrodynamics/OrbitDetermination/EstimatableParameters/rotationalPrecessionRate.h"
#include "Astrodynamics/OrbitDetermination/ObservationPartials/createPositionPartials.h"
#include "Astrodynamics/OrbitDetermination/ObservationPartials/numericalObservationPartial.h"
#include "SimulationSetup/createBodyDeformationModel.h"
#include "SimulationSetup/createGroundStations.h"
#include "SimulationSetup/createBodies.h"
#include "SimulationSetup/defaultBodies.h"
#include "SimulationSetup/createEstimatableParameters.h"


namespace tudat
{
namespace unit_tests
{

using namespace tudat::bodies;
using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::site_displacements;
using namespace tudat::spice_interface;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::observation_partials;

BOOST_AUTO_TEST_SUITE( test_position_partials)

BOOST_AUTO_TEST_CASE( testPositionPartials )
{

    //Load spice kernels.
    std::string kernelsPath = input_output::getDataFilesRootPath( ) + "SpiceKernels/";
    loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    loadSpiceKernelInTudat( kernelsPath + "de421.bsp");
    loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");
    loadSpiceKernelInTudat( kernelsPath + "mar097.bsp");

    // Specify initial time
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = 1.2E7;
    double buffer = 1000.0;

    // Create bodies.
    NamedBodyMap bodyMap;
    bodyMap[ "Earth" ] = boost::make_shared< CelestialBody >( );
    bodyMap[ "Moon" ] = boost::make_shared< CelestialBody >( );
    bodyMap[ "Sun" ] = boost::make_shared< CelestialBody >( );
    bodyMap[ "Mars" ] = boost::make_shared< CelestialBody >( );
    bodyMap[ "Phobos" ] = boost::make_shared< CelestialBody >( );

    boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Earth" ] )->setShapeModel( boost::make_shared< SphericalBodyShapeModel >(
                                                                                           spice_interface::getAverageRadius( "Earth" ) ) );
    boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Mars" ] )->setShapeModel( boost::make_shared< SphericalBodyShapeModel >(
                                                                                           spice_interface::getAverageRadius( "Mars" ) ) );
    boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Phobos" ] )->setShapeModel( boost::make_shared< SphericalBodyShapeModel >(
                                                                                            spice_interface::getAverageRadius( "Phobos" ) ) );

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
    bodyMap[ "Phobos" ]->setEphemeris(
                boost::make_shared< ConstantEphemeris >(
                    getBodyCartesianStateAtEpoch( "Phobos", "Mars", "ECLIPJ2000", "NONE", 0.0 ),
                    "SSB", "ECLIPJ2000" ) );

    boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Sun" ] )->setGravityFieldModel(
                boost::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Sun" ) ) );
    boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Moon" ] )->setGravityFieldModel(
                boost::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Moon" ) ) );
    boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Earth" ] )->setGravityFieldModel(
                boost::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Earth" ) ) );

    std::vector< Eigen::Vector3d > dummyMeanDeformingBodyPositionVectors;
    dummyMeanDeformingBodyPositionVectors.push_back(
                ( Eigen::Vector3d( )<<400E6, 0.0, 0.0 ).finished( ) );
    dummyMeanDeformingBodyPositionVectors.push_back(
                ( Eigen::Vector3d( )<<0.0, 1.5E11, 0.0 ).finished( ) );

    boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Earth" ] )->setBasicTidalBodyDeformation(
                boost::dynamic_pointer_cast< BasicTidalBodyDeformation >(
                    createBodyDeformationModel(
                        getDefaultBodyDeformationSettings(
                            "Earth", initialEphemerisTime, finalEphemerisTime )[ 0 ], "Earth", bodyMap ) ) );
    boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Earth" ] )->getBasicTidalBodyDeformation( )->setMeanVectorToPerturbingBodies(
                dummyMeanDeformingBodyPositionVectors );

    boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Earth" ] )->setRotationalEphemeris(
                createRotationModel( getDefaultRotationModelSettings(
                                         "Earth", initialEphemerisTime, finalEphemerisTime ), "Earth", bodyMap ) );
    boost::shared_ptr< SimpleRotationalEphemeris > simpleEarthRotationModel =
            boost::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                bodyMap[ "Earth" ]->getRotationalEphemeris( ) );

    bodyMap[ "Mars" ]->setRotationalEphemeris( createRotationModel( getDefaultRotationModelSettings( "Mars", initialEphemerisTime,
                                                                                                     finalEphemerisTime, full ), "Mars",  bodyMap ) );
    bodyMap[ "Phobos" ]->setRotationalEphemeris( createRotationModel( getDefaultRotationModelSettings( "Phobos", initialEphemerisTime,
                                                                                                       finalEphemerisTime, full ), "Phobos", bodyMap ) );

    boost::shared_ptr< DirectlyPerturbedRotationModel > phobosDirectRotationModel = boost::dynamic_pointer_cast< DirectlyPerturbedRotationModel >(
                bodyMap[ "Phobos" ]->getRotationalEphemeris( ) );
    boost::shared_ptr< PlanetaryOrientationAngleCalculator > marsRotationAngleCalculator = bodyMap[ "Mars" ]->getPlanetaryOrientationAnglesCalculator( );
    boost::shared_ptr< PlanetaryRotationModel > marsRotationModel = boost::dynamic_pointer_cast< PlanetaryRotationModel >(
                bodyMap[ "Mars" ]->getRotationalEphemeris( ) );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );




    // Define and create ground stations.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.push_back( std::make_pair( "Earth", "Graz" ) );



    LinkEnds linkEnds;
    linkEnds[ observed_body ] = groundStations[ 0 ];
    createGroundStations( bodyMap, groundStations );


    Eigen::Vector3d phobosLanderGeodeticCoordinates = ( Eigen::Vector3d( )<< 0.0, -.25, -.3 ).finished( );
    boost::shared_ptr< NominalGroundStationState > stationState =
            boost::make_shared< GeodeticNominalGroundStationState >(
                phobosLanderGeodeticCoordinates.y( ),
                phobosLanderGeodeticCoordinates.z( ),
                phobosLanderGeodeticCoordinates.x( ),
                boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Phobos" ] )->getShapeModel( ),
                Eigen::Vector3d::Zero( ),
                basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    boost::shared_ptr< PointingAnglesCalculator > pointingAnglesCalculator =
            boost::make_shared< PointingAnglesCalculator >(
                boost::bind( &RotationalEphemeris::getRotationToTargetFrame, bodyMap.at( "Phobos" )->getRotationalEphemeris( ), _1 ),
                boost::bind( &NominalGroundStationState::getRotationFromBodyFixedToTopocentricFrame, stationState, _1 ) );

    boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Phobos" ] )->addGroundStation(
                "Lander", boost::make_shared< GroundStation >(
                    stationState, pointingAnglesCalculator,
                    "Lander" ) );

    groundStations.push_back( std::make_pair( "Phobos", "Lander" ) );
    LinkEnds linkEnds2;
    linkEnds2[ observed_body ] = groundStations[ 1 ];


    Eigen::Vector3d marsLanderGeodeticCoordinates = ( Eigen::Vector3d( )<< 0.0, .35, -.214 ).finished( );
    boost::shared_ptr< NominalGroundStationState > marsStationState =
            boost::make_shared< GeodeticNominalGroundStationState >(
                marsLanderGeodeticCoordinates.y( ),
                marsLanderGeodeticCoordinates.z( ),
                marsLanderGeodeticCoordinates.x( ),
                boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Mars" ] )->getShapeModel( ),
                Eigen::Vector3d::Zero( ),
                basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    boost::shared_ptr< PointingAnglesCalculator > marsPointingAnglesCalculator =
            boost::make_shared< PointingAnglesCalculator >(
                boost::bind( &RotationalEphemeris::getRotationToTargetFrame, bodyMap.at( "Mars" )->getRotationalEphemeris( ), _1 ),
                boost::bind( &NominalGroundStationState::getRotationFromBodyFixedToTopocentricFrame, marsStationState, _1 ) );

    boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Mars" ] )->addGroundStation(
                "Lander2", boost::make_shared< GroundStation >(
                    stationState, marsPointingAnglesCalculator,
                    "Lander2" ) );

    groundStations.push_back( std::make_pair( "Mars", "Lander2" ) );
    LinkEnds linkEnds3;
    linkEnds3[ observed_body ] = groundStations[ 2 ];

    boost::shared_ptr< GroundStation > receivingGroundStation =
            boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Earth" ] )->getGroundStation( "Graz" );

    setSingleBodyGroundStationPositionVariationFunctions( boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Earth" ] ), basic_solid_body,
                                                          initialEphemerisTime, finalEphemerisTime );

    boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Earth" ] )->updateGroundStationStateHistory( );

    //Create parameter objects.
    std::string earthName = "Earth";
    boost::shared_ptr< H2DisplacementLoveNumber > earthH2LoveNumber = boost::make_shared< H2DisplacementLoveNumber  >(
                boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Earth" ] )->getBasicTidalBodyDeformation( ), "Earth" );
    boost::shared_ptr< L2DisplacementShidaNumber > earthL2LoveNumber = boost::make_shared< L2DisplacementShidaNumber  >(
                boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Earth" ] )->getBasicTidalBodyDeformation( ), "Earth" );
    boost::shared_ptr< RotationRate > earthRotationRate = boost::make_shared< RotationRate >(
                boost::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                    bodyMap[ "Earth" ]->getRotationalEphemeris( ) ), earthName);
    boost::shared_ptr< ConstantRotationalOrientation > earthPolePosition = boost::make_shared< ConstantRotationalOrientation >(
                boost::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                    bodyMap[ "Earth" ]->getRotationalEphemeris( ) ), earthName );
    boost::shared_ptr< GroundStationPosition > grazPosition = boost::make_shared< GroundStationPosition >(
                boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Earth" ] )->getGroundStation( "Graz" )->getNominalStationState( ),
                "Earth", "Graz" );

    boost::shared_ptr< RotationModelPeriodicVariationAmplitudes > rightAscensionAmplitudes =
            boost::make_shared< RotationModelPeriodicVariationAmplitudes >(
                phobosDirectRotationModel, phobosDirectRotationModel->getRightAscensionTermPeriods( ), right_ascension_angle, "Phobos" );
    boost::shared_ptr< RotationModelPeriodicVariationAmplitudes > declinationAmplitudes =
            boost::make_shared< RotationModelPeriodicVariationAmplitudes >(
                phobosDirectRotationModel, phobosDirectRotationModel->getDeclinationTermPeriods( ), declination_angle, "Phobos" );
    boost::shared_ptr< RotationModelPeriodicVariationAmplitudes > primeMeridianAmplitudes =
            boost::make_shared< RotationModelPeriodicVariationAmplitudes >(
                phobosDirectRotationModel, phobosDirectRotationModel->getPrimeMeridianTermPeriods( ), prime_meridian_angle, "Phobos" );

    boost::shared_ptr< RotationModelPolynomialVariations > rightAscensionPolynomialAmplitudes =
            boost::make_shared< RotationModelPolynomialVariations >(
                phobosDirectRotationModel, boost::assign::list_of( 0 )( 1 ), right_ascension_angle, "Phobos" );
    boost::shared_ptr< RotationModelPolynomialVariations > declinationPolynomialAmplitudes =
            boost::make_shared< RotationModelPolynomialVariations >(
                phobosDirectRotationModel, boost::assign::list_of( 0 )( 1 ), declination_angle, "Phobos" );
    boost::shared_ptr< RotationModelPolynomialVariations > primeMeridianPolynomialAmplitudes =
            boost::make_shared< RotationModelPolynomialVariations >(
                phobosDirectRotationModel, boost::assign::list_of( 0 )( 1 )( 2 ), prime_meridian_angle, "Phobos" );

    boost::shared_ptr< RotationalPrecessionRate > rotationalPrecessionRate =
            boost::make_shared< RotationalPrecessionRate >( marsRotationAngleCalculator, "Mars" );


    // Create explicit range partial objects.
    boost::shared_ptr< PositionPartial > partialObjectWrtReceiverPosition =
            createPositionPartialsWrtBodyPosition( linkEnds, bodyMap, "Earth" ).begin( )->second;

    boost::shared_ptr< PositionPartial > partialObjectWrtReceiverRotationRate =
            createPositionPartialsWrtParameter(
                linkEnds, bodyMap, earthRotationRate ).begin( )->second;
    boost::shared_ptr< PositionPartial > partialObjectWrtReceiverH2 =
            createPositionPartialsWrtParameter(
                linkEnds, bodyMap, earthH2LoveNumber ).begin( )->second;
    boost::shared_ptr< PositionPartial > partialObjectWrtReceiverL2 =
            createPositionPartialsWrtParameter(
                linkEnds, bodyMap, earthL2LoveNumber ).begin( )->second;
    boost::shared_ptr< PositionPartial > partialObjectWrtReceiverPolePosition =
            createPositionPartialsWrtParameter(
                linkEnds, bodyMap, earthPolePosition ).begin( )->second;
    boost::shared_ptr< PositionPartial > partialObjectWrtReceiverGroundStationPosition =
            createPositionPartialsWrtParameter(
                linkEnds, bodyMap, grazPosition ).begin( )->second;

    boost::shared_ptr< PositionPartial > partialObjectWrtRightAscensionAmplitudes =
            createPositionPartialsWrtParameter(
                linkEnds2, bodyMap, rightAscensionAmplitudes ).begin( )->second;
    boost::shared_ptr< PositionPartial > partialObjectWrtDeclinationAmplitudes =
            createPositionPartialsWrtParameter(
                linkEnds2, bodyMap, declinationAmplitudes ).begin( )->second;
    boost::shared_ptr< PositionPartial > partialObjectWrtPrimeMeridianAmplitudes =
            createPositionPartialsWrtParameter(
                linkEnds2, bodyMap, primeMeridianAmplitudes ).begin( )->second;

    boost::shared_ptr< PositionPartial > partialObjectWrtRightAscensionPolynomialAmplitudes =
            createPositionPartialsWrtParameter(
                linkEnds2, bodyMap, rightAscensionPolynomialAmplitudes ).begin( )->second;
    boost::shared_ptr< PositionPartial > partialObjectWrtDeclinationPolynomialAmplitudes =
            createPositionPartialsWrtParameter(
                linkEnds2, bodyMap, declinationPolynomialAmplitudes ).begin( )->second;
    boost::shared_ptr< PositionPartial > partialObjectWrtPrimeMeridianPolynomialAmplitudes =
            createPositionPartialsWrtParameter(
                linkEnds2, bodyMap, primeMeridianPolynomialAmplitudes ).begin( )->second;
    boost::shared_ptr< PositionPartial > partialObjectWrtMarsPrecessionRate =
            createPositionPartialsWrtParameter(
                linkEnds3, bodyMap, rotationalPrecessionRate ).begin( )->second;

    // Calculate transmission/reception times and states
    basic_mathematics::Vector6d currentState;
    double receptionTime = 1.1E7;
    currentState = bodyMap.at( "Earth" )->getStateInBaseFrameFromEphemeris( receptionTime );

    double currentTime = receptionTime;

    Eigen::MatrixXd partialWrtReceiverPosition =
            partialObjectWrtReceiverPosition->calculatePartial( currentState, currentTime );

    Eigen::MatrixXd partialWrtReceiverRotationrate =
            partialObjectWrtReceiverRotationRate->calculatePartial( currentState, currentTime );

    Eigen::MatrixXd partialWrtReceiverL2 =
            partialObjectWrtReceiverL2->calculatePartial( currentState, currentTime );

    Eigen::MatrixXd partialWrtReceiverH2 =
            partialObjectWrtReceiverH2->calculatePartial( currentState, currentTime );

    Eigen::MatrixXd partialWrtReceiverPolePosition =
            partialObjectWrtReceiverPolePosition->calculatePartial( currentState, currentTime );

    Eigen::MatrixXd partialWrtGroundStationPosition =
            partialObjectWrtReceiverGroundStationPosition->calculatePartial( currentState, currentTime );

    Eigen::MatrixXd partialWrtRightAscensionAmplitudes =
            partialObjectWrtRightAscensionAmplitudes->calculatePartial( currentState, currentTime );

    Eigen::MatrixXd partialWrtDeclinationAmplitudes =
            partialObjectWrtDeclinationAmplitudes->calculatePartial( currentState, currentTime );

    Eigen::MatrixXd partialWrtPrimeMeridianAmplitudes =
            partialObjectWrtPrimeMeridianAmplitudes->calculatePartial( currentState, currentTime );

    Eigen::MatrixXd partialWrtRightAscensionPolynomialAmplitudes =
            partialObjectWrtRightAscensionPolynomialAmplitudes->calculatePartial( currentState, currentTime );

    Eigen::MatrixXd partialWrtDeclinationPolynomialAmplitudes =
            partialObjectWrtDeclinationPolynomialAmplitudes->calculatePartial( currentState, currentTime );

    Eigen::MatrixXd partialWrtPrimeMeridianPolynomialAmplitudes =
            partialObjectWrtPrimeMeridianPolynomialAmplitudes->calculatePartial( currentState, currentTime );

    Eigen::MatrixXd partialWrtMarsPrecessionRate =
            partialObjectWrtMarsPrecessionRate->calculatePartial( currentState, currentTime );


    boost::function< Eigen::VectorXd( const double ) > observationFunctionAtReception =
            boost::bind( &Ephemeris::getPosition, createReferencePointEphemeris< double, double >(
                             bodyMap.at( "Earth" ), bodyMap.at( "Earth" )->getRotationalEphemeris( ),
                             boost::bind( &GroundStation::getStateInPlanetFixedFrame< double >,
                                          boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Earth" ] )->getGroundStation( "Graz" ), _1 ) ), _1,
                         basic_astrodynamics::JULIAN_DAY_ON_J2000 );


    double rotationRateVariation = 1.0E-10;
    double loveNumberVariation = 1.0E-5;

    Eigen::Vector3d bodyPositionVariation;
    bodyPositionVariation << 10.0, 10.0, 10.0;

    // Calculate numerical partials w.r.t. Earth state.
    boost::shared_ptr< ConstantEphemeris > earthEphemeris = boost::dynamic_pointer_cast< ConstantEphemeris >(
                bodyMap[ "Earth" ]->getEphemeris( ) );
    basic_mathematics::Vector6d earthUnperturbedState = earthEphemeris->getCartesianStateFromEphemeris( 0.0 );
    basic_mathematics::Vector6d perturbedEarthState;
    Eigen::Matrix< double, 3, 3 > numericalPartialWrtReceiverPosition = Eigen::Matrix< double, 3, 3 >::Zero( );

    for( int i = 0; i < 3; i++ )
    {
        perturbedEarthState = earthUnperturbedState;
        perturbedEarthState( i ) += bodyPositionVariation( i );
        earthEphemeris->updateConstantState( perturbedEarthState );
        Eigen::Vector3d upPerturbedPosition = observationFunctionAtReception( currentTime );

        perturbedEarthState = earthUnperturbedState;
        perturbedEarthState( i ) -= bodyPositionVariation( i );
        earthEphemeris->updateConstantState( perturbedEarthState );
        Eigen::Vector3d downPerturbedPosition = observationFunctionAtReception( currentTime );

        numericalPartialWrtReceiverPosition.block( 0, i, 3, 1 ) = ( upPerturbedPosition - downPerturbedPosition ) /
                ( 2.0 * bodyPositionVariation( i ) );
    }
    earthEphemeris->updateConstantState( earthUnperturbedState );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtReceiverPosition, numericalPartialWrtReceiverPosition, 1.0E-4 );



    boost::function< void( ) > earthDisplacementUpdateFunction = boost::bind(
                &CelestialBody::updateGroundStationStateHistory, boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Earth" ] ) );

    Eigen::Vector3d numericalPartialWrtReceiverH2LoveNumber = calculateNumericalObservationParameterPartial(
                earthH2LoveNumber, 1.0, observationFunctionAtReception, receptionTime,
                earthDisplacementUpdateFunction );
    Eigen::Vector3d numericalPartialWrtReceiverL2LoveNumber =calculateNumericalObservationParameterPartial(
                earthL2LoveNumber, 1.0, observationFunctionAtReception, receptionTime,
                earthDisplacementUpdateFunction );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtReceiverH2, numericalPartialWrtReceiverH2LoveNumber, 2.0E-4 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtReceiverL2, numericalPartialWrtReceiverL2LoveNumber, 2.0E-4 );

    Eigen::Vector3d numericalPartialWrtReceiverRotationrate = calculateNumericalObservationParameterPartial(
                earthRotationRate, 1.0E-8, observationFunctionAtReception,
                receptionTime );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtReceiverRotationrate, numericalPartialWrtReceiverRotationrate, 1.0E-4 );



    Eigen::VectorXd polePositionPerturbation = ( Eigen::Vector2d( )<<1.0E-4, 1.0E-4 ).finished( );
    Eigen::MatrixXd numericalPartialWrtReceiverPolePosition = calculateNumericalObservationParameterPartial(
                earthPolePosition, polePositionPerturbation, observationFunctionAtReception, receptionTime );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtReceiverPolePosition, numericalPartialWrtReceiverPolePosition, 1.0E-4 );


    Eigen::VectorXd groundStationPerturbation = ( Eigen::Vector3d( )<<100.0, 100.0, 100.0 ).finished( );
    Eigen::MatrixXd numericalPartialWrtReceiverGroundStationPosition = calculateNumericalObservationParameterPartial(
                grazPosition, groundStationPerturbation, observationFunctionAtReception, receptionTime );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtGroundStationPosition, numericalPartialWrtReceiverGroundStationPosition, 1.0E-4 );


    // Calculate partials of Phobos laser position
    boost::function< Eigen::VectorXd( const double ) > phobosObservationFunctionAtReception =
            boost::bind( &Ephemeris::getPosition, createReferencePointEphemeris< double, double >(
                             bodyMap.at( "Phobos" ), bodyMap.at( "Phobos" )->getRotationalEphemeris( ),
                             boost::bind( &GroundStation::getStateInPlanetFixedFrame< double >,
                                          boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Phobos" ] )->getGroundStation( "Lander" ), _1 ) ), _1,
                         basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    boost::function< void( ) > rotationModelUpdateFunction = boost::bind( &DirectlyPerturbedRotationModel::update, phobosDirectRotationModel,
                                                                          receptionTime );

    Eigen::VectorXd rightAscensionAmplitudePerturbations = rightAscensionAmplitudes->getParameterValue( );
    for( unsigned int i = 0; i < rightAscensionAmplitudePerturbations.rows( ); i++ )
    {
        rightAscensionAmplitudePerturbations( i ) = 1.0E-3 * mathematical_constants::PI / 180.0;
    }
    Eigen::MatrixXd numericalPartialWrtRightAscensionAmplitudes = calculateNumericalObservationParameterPartial(
                rightAscensionAmplitudes, rightAscensionAmplitudePerturbations, phobosObservationFunctionAtReception, receptionTime,
                rotationModelUpdateFunction );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtRightAscensionAmplitudes, numericalPartialWrtRightAscensionAmplitudes, 1.0E-4 );


    Eigen::VectorXd declinationAmplitudePerturbations = declinationAmplitudes->getParameterValue( );
    for( unsigned int i = 0; i < declinationAmplitudePerturbations.rows( ); i++ )
    {
        declinationAmplitudePerturbations( i ) = 1.0E-3 * mathematical_constants::PI / 180.0;
    }
    Eigen::MatrixXd numericalPartialWrtDeclinationAmplitudes = calculateNumericalObservationParameterPartial(
                declinationAmplitudes, declinationAmplitudePerturbations, phobosObservationFunctionAtReception, receptionTime,
                rotationModelUpdateFunction );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtDeclinationAmplitudes, numericalPartialWrtDeclinationAmplitudes, 1.0E-4 );

    Eigen::VectorXd primeMeridianAmplitudePerturbations = primeMeridianAmplitudes->getParameterValue( );
    for( unsigned int i = 0; i < primeMeridianAmplitudePerturbations.rows( ); i++ )
    {
        primeMeridianAmplitudePerturbations( i ) = 1.0E-3 * mathematical_constants::PI / 180.0;
    }
    Eigen::MatrixXd numericalPartialWrtPrimeMeridianAmplitudes = calculateNumericalObservationParameterPartial(
                primeMeridianAmplitudes, primeMeridianAmplitudePerturbations, phobosObservationFunctionAtReception, receptionTime,
                rotationModelUpdateFunction );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPrimeMeridianAmplitudes, numericalPartialWrtPrimeMeridianAmplitudes, 1.0E-4 );


    // Calculate partials of Mars lander position
    boost::function< Eigen::VectorXd( const double ) > marsObservationFunctionAtReception =
            boost::bind( &Ephemeris::getPosition, createReferencePointEphemeris< double, double >(
                             bodyMap.at( "Mars" ), bodyMap.at( "Mars" )->getRotationalEphemeris( ),
                             boost::bind( &GroundStation::getStateInPlanetFixedFrame< double >,
                                          boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Mars" ] )->getGroundStation( "Lander2" ), _1 ) ), _1,
                         basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    //boost::function< void( ) > rotationModelUpdateFunction = boost::bind( &PlanetaryOrientationAngleCalculator::update, phobosDirectRotationModel,
    //                                                                      receptionTime );

    Eigen::MatrixXd numericalPartialWrtMarsPrecessionRate = calculateNumericalObservationParameterPartial(
                rotationalPrecessionRate, 100.0 * rotationalPrecessionRate->getParameterValue( ), marsObservationFunctionAtReception,
                receptionTime );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtMarsPrecessionRate, numericalPartialWrtMarsPrecessionRate, 1.0E-4 );
}

BOOST_AUTO_TEST_CASE( testTidallyLockedLibrationAmplitudePartial )
{
    using namespace simulation_setup;
    using namespace input_output;

    std::string kernelsPath = input_output::getDataFilesRootPath( ) + "SpiceKernels/";
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "jup310.bsp");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00010.tpc");

    std::map< std::string, BodySettingLevel > bodyNames;

    bodyNames[ "Sun" ] = simple;
    bodyNames[ "Jupiter" ] = simple;
    bodyNames[ "Io" ] = full;

    // Specify initial time
    double initialTime = 1.0E7;
    double finalTime = 1.05E7;

    // Test without librations
    {
        // Create bodies needed in simulation
        std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
                getDefaultBodySettings< double, double >( bodyNames, initialTime, finalTime );

        simulation_setup::NamedBodyMap bodyMap = createCelestialBodies( bodySettings );

        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );
        boost::shared_ptr< RotationalEphemeris > rotationModel = bodyMap[ "Io" ]->getRotationalEphemeris( );

        std::map< std::string, std::pair< int, int > > librationType;
        librationType[ "Io" ] = std::make_pair( 5, 1 );
        boost::shared_ptr< estimatable_parameters::EstimatableParameterSettings > parametersToEstimate =
                boost::make_shared< TidallyLockedLibrationAmplitudeSettings >( librationType, "Io" );

        boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > librationAmplitude =
                createDoubleParameterToEstimate( parametersToEstimate, bodyMap );
        boost::shared_ptr< RotationMatrixPartial > rotationMatrixPartial = createRotationMatrixPartialsWrtParameter(
                    bodyMap, librationAmplitude );

        double librationPerturbation = 1.0E-8;
        Eigen::Matrix3d upPerturbedRotation, downPerturbedRotation;

        double nominalLibration = librationAmplitude->getParameterValue( );

        librationAmplitude->setParameterValue( nominalLibration + librationPerturbation );
        upPerturbedRotation =  ( rotationModel->getRotationToBaseFrame( 1.025E7 ) );

        librationAmplitude->setParameterValue( nominalLibration - librationPerturbation );
        downPerturbedRotation = ( rotationModel->getRotationToBaseFrame( 1.025E7 ) );

        librationAmplitude->setParameterValue( nominalLibration );

        Eigen::Matrix3d numericalPartial = ( upPerturbedRotation - downPerturbedRotation ) / ( 2.0 * librationPerturbation );
        Eigen::Matrix3d analyticalPartial = rotationMatrixPartial->calculatePartialOfRotationMatrixWrParameter( 1.025E7 ).at( 0 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( numericalPartial.block( 0, 0, 3, 1 ) ), ( analyticalPartial.block( 0, 0, 3, 1 ) ), 1.0E-7 );
        BOOST_CHECK_SMALL( numericalPartial( 0, 2 ), 1.0E-10 );
        BOOST_CHECK_SMALL( numericalPartial( 1, 2 ), 1.0E-10 );
        BOOST_CHECK_SMALL( numericalPartial( 2, 2 ), 1.0E-10 );

        std::cout<<( numericalPartial - analyticalPartial ).cwiseQuotient( analyticalPartial)<<std::endl;


    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat





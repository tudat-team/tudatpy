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


#include <limits>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/astro/ephemerides/keplerEphemeris.h"
#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"
#include "tudat/simulation/estimation.h"

namespace tudat
{
namespace unit_tests
{

//Using declarations.
using namespace tudat::ephemerides;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::propagators;
using namespace tudat::coordinate_conversions;
using namespace tudat::reference_frames;
using namespace tudat::observation_models;
using namespace tudat::estimatable_parameters;
using namespace tudat::orbit_determination;


BOOST_AUTO_TEST_SUITE( test_rotational_dynamics_estimation )

//! Create system of bodies to be used in estimation test
SystemOfBodies getTestBodyMap( const double phobosSemiMajorAxis,
                             const bool useSymmetricEquator = 0 )
{
    SystemOfBodies bodies = SystemOfBodies( "Mars", "ECLIPJ2000" );

    // Create Mars object
    bodies.createEmptyBody( "Mars", false );
    bodies.at( "Mars" )->setEphemeris( createBodyEphemeris(
                                         getDefaultEphemerisSettings( "Mars" ), "Mars" ) );
    std::shared_ptr< SphericalHarmonicsGravityFieldSettings > marsGravityFieldSettings =
            std::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >(
                getDefaultGravityFieldSettings( "Mars", TUDAT_NAN, TUDAT_NAN ) );
    marsGravityFieldSettings->resetGravitationalParameter( spice_interface::getBodyGravitationalParameter( "Mars" ) );
    bodies.at( "Mars" )->setGravityFieldModel(
                createGravityFieldModel( marsGravityFieldSettings, "Mars", bodies ) );
    bodies.at( "Mars" )->setRotationalEphemeris(
                createRotationModel( getDefaultRotationModelSettings( "Mars", TUDAT_NAN, TUDAT_NAN ), "Mars" ) );

    // Create Mars object
    bodies.createEmptyBody( "Earth" );
    bodies.at( "Earth" )->setEphemeris( createBodyEphemeris(
                                          getDefaultEphemerisSettings( "Earth" ), "Earth" ) );

    // Create Phobos object
    bodies.createEmptyBody( "Phobos" );

    // Set Phobos inertia
    Eigen::Matrix3d phobosInertiaTensor = Eigen::Matrix3d::Zero( );
    phobosInertiaTensor( 0, 0 ) = 0.3615;
    phobosInertiaTensor( 1, 1 ) = 0.4265;
    phobosInertiaTensor( 2, 2 ) = 0.5024;
    if( useSymmetricEquator )
    {
        phobosInertiaTensor( 0, 0 ) = phobosInertiaTensor( 1, 1 );
    }
    double phobosReferenceRadius = 11.27E3;
    double phobosMass = 1.0659E16;
    phobosInertiaTensor *= (phobosReferenceRadius * phobosReferenceRadius * phobosMass );
    bodies.at( "Phobos" )->setBodyInertiaTensor( phobosInertiaTensor );

    // Set Phobos shape
    bodies.at( "Phobos" )->setShapeModel(
                std::make_shared< SphericalBodyShapeModel >( 15.0E3 ) );

    // Compute and set Phobos gravity field
    double phobosGravitationalParameter = phobosMass * physical_constants::GRAVITATIONAL_CONSTANT;
    Eigen::MatrixXd phobosCosineGravityFieldCoefficients = Eigen::Matrix3d::Zero( ),
            phobosSineGravityFieldCoefficients = Eigen::Matrix3d::Zero( );
    double phobosScaledMeanMomentOfInertia;
    gravitation::getDegreeTwoSphericalHarmonicCoefficients(
                phobosInertiaTensor, phobosGravitationalParameter, phobosReferenceRadius, true,
                phobosCosineGravityFieldCoefficients, phobosSineGravityFieldCoefficients, phobosScaledMeanMomentOfInertia );
    bodies.at( "Phobos" )->setGravityFieldModel(
                std::make_shared< gravitation::SphericalHarmonicsGravityField >(
                    phobosGravitationalParameter, phobosReferenceRadius, phobosCosineGravityFieldCoefficients,
                    phobosSineGravityFieldCoefficients, "Phobos_Fixed",
                    std::bind( &Body::setBodyInertiaTensorFromGravityFieldAndExistingMeanMoment, bodies.at( "Phobos" ), true ) ) );


    // Set Phobos dummy rotational ephemeris
    Eigen::Quaterniond noRotationQuaternion = Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );
    Eigen::Matrix< double, 7, 1 > unitRotationState = Eigen::Matrix< double, 7, 1 >::Zero( );
    unitRotationState( 0 ) = noRotationQuaternion.w( );
    unitRotationState( 1 ) = noRotationQuaternion.x( );
    unitRotationState( 2 ) = noRotationQuaternion.y( );
    unitRotationState( 3 ) = noRotationQuaternion.z( );

    std::map< double, Eigen::Matrix< double, 7, 1 > > dummyRotationMap;
    dummyRotationMap[ -1.0E100 ] = unitRotationState;
    dummyRotationMap[ 1.0E100 ] = unitRotationState;

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 7, 1 > > > dummyInterpolator =
            std::make_shared< interpolators::LinearInterpolator< double, Eigen::Matrix< double, 7, 1 > > >( dummyRotationMap );
    bodies.at( "Phobos" )->setRotationalEphemeris( std::make_shared< TabulatedRotationalEphemeris< double, double > >(
                                                     dummyInterpolator, "ECLIPJ2000", "Phobos_Fixed" ) );

    // Set Phobos on Kepler orbit
    Eigen::Vector6d phobosKeplerElements = Eigen::Vector6d::Zero( );
    phobosKeplerElements( 0 ) = phobosSemiMajorAxis;
    bodies.at( "Phobos" )->setEphemeris(
                ephemerides::getTabulatedEphemeris(
                    std::make_shared< ephemerides::KeplerEphemeris >(
                        phobosKeplerElements, 0.0, spice_interface::getBodyGravitationalParameter( "Mars" ),
                        "Mars", "ECLIPJ2000" ), -3600.0, 120.0 * 86400.0 + 3600.0, 120.0,
                    std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) ) );

    return bodies;
}

//! Test if Phobos rotational dynamics is correctly estimation from lander tracking data
BOOST_AUTO_TEST_CASE( test_RotationalDynamicsEstimationFromLanderData )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Retrieve list of body objects.
    SystemOfBodies bodies = getTestBodyMap( 9376.0E3, 0 );
    createGroundStation( bodies.at( "Phobos" ), "Lander", ( Eigen::Vector3d( ) << 0.1, 0.35, 0.0 ).finished( ), geodetic_position );

    // Define time range of test.
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 60.0 * 86400.0;

    // Define mean motion (equal to rotation rate).
    double phobosSemiMajorAxis = 9376.0E3;
    double meanMotion = std::sqrt( getBodyGravitationalParameter( "Mars" ) /
                                   std::pow( phobosSemiMajorAxis, 3.0 ) );

    // Define initial rotational state (slight-perturbation from tidally locked)
    double initialAnglePerturbation = 1.0E-3;
    double initialRotationRatePerturbation = 1.0E-2;
    Eigen::Quaterniond nominalInitialRotation =
            Eigen::Quaterniond( Eigen::AngleAxisd( -initialAnglePerturbation, Eigen::Vector3d::UnitX( ) ) );
    Eigen::VectorXd systemInitialState = Eigen::VectorXd::Zero( 7 );
    systemInitialState.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat( nominalInitialRotation );
    systemInitialState( 6 ) = meanMotion * ( 1.0 + initialRotationRatePerturbation );
//    Eigen::Matrix3d phobosInertiaTensor = bodies.at( "Phobos" )->getBodyInertiaTensor( );

    // Create torque models
    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "Phobos" );
    SelectedTorqueMap torqueMap;
    torqueMap[ "Phobos" ][ "Mars" ].push_back( std::make_shared< TorqueSettings >( second_order_gravitational_torque ) );
    basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap(
                bodies, torqueMap, bodiesToIntegrate );

    // Define integrator settings.
    double timeStep = 240.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, initialEphemerisTime, timeStep );

    // Define propagator settings.
    std::shared_ptr< RotationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< RotationalStatePropagatorSettings< double > >
            ( torqueModelMap, bodiesToIntegrate, systemInitialState, std::make_shared< PropagationTimeTerminationSettings >(
                  finalEphemerisTime ) );

    // Define link ends and observables
    std::vector< LinkEnds > linkEndsList;
    LinkEnds currentLinkEnds;
    currentLinkEnds[ transmitter ] = std::make_pair( "Earth", "" );
    currentLinkEnds[ receiver ] = std::make_pair( "Phobos", "Lander" );
    linkEndsList.push_back( currentLinkEnds );
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    //linkEndsPerObservable[ one_way_range ].push_back( linkEndsList[ 0 ] );

    LinkEnds currentLinkEnds2;
    currentLinkEnds2[ observed_body ] = std::make_pair( "Phobos", "" );
    linkEndsPerObservable[ euler_angle_313_observable ].push_back( currentLinkEnds2 );

    // Create parameters to estimate
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialStateParameterSettings< double >( propagatorSettings, bodies );

    std::vector< std::pair< int, int > > blockIndices;
    blockIndices.push_back( std::make_pair( 2, 2 ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  blockIndices, "Phobos", spherical_harmonics_cosine_coefficient_block ) );
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodies );
    printEstimatableParameterEntries( parametersToEstimate );

    // Create observation settings
    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;

        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                                                               currentObservable, currentLinkEndsList.at( i ) ) );
        }
    }

    // Create orbit determination object
    OrbitDeterminationManager< double, double > orbitDeterminationManager =
            OrbitDeterminationManager< double, double >(
                bodies, parametersToEstimate, observationSettingsList,
                integratorSettings, propagatorSettings );

    // Deifne observation times
    std::vector< double > observationTimes;
    double currentTime = initialEphemerisTime + 1800.0;
    double observationTimeStep = 60.0;
    while( currentTime < finalEphemerisTime - 1800.0 )
    {
        observationTimes.push_back( currentTime );
        currentTime += observationTimeStep;
    }

    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;

    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            measurementSimulationInput.push_back(
                        std::make_shared< TabulatedObservationSimulationSettings< > >(
                            currentObservable, currentLinkEndsList.at( i ), observationTimes, observed_body ) );
        }
    }

    // Simulate observations
    std::shared_ptr< ObservationCollection< > > observationsAndTimes = simulateObservations< double, double >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    // Perturb parameter estimate
    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< double >( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    int numberOfParameters = initialParameterEstimate.rows( );
    initialParameterEstimate( 1 ) -= 1.0E-7;
    initialParameterEstimate( 2 ) -= 1.0E-7;
    initialParameterEstimate( 3 ) += 1.0E-7;

    initialParameterEstimate.segment( 0, 4 ).normalize( );
    initialParameterEstimate( 5 ) += 1.0E-9;
    initialParameterEstimate( 7 ) += 1.0E-12;

    // Define estimation input
    std::shared_ptr< PodInput< double, double  > > podInput =
            std::make_shared< PodInput< double, double > >(
                observationsAndTimes, numberOfParameters,
                Eigen::MatrixXd::Zero( numberOfParameters, numberOfParameters ), initialParameterEstimate - truthParameters );

    // Perform estimation
    std::shared_ptr< PodOutput< double > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, std::make_shared< EstimationConvergenceChecker >( 6 ) );

    // Check residual size (sub-mm over >1 AU)
    BOOST_CHECK_SMALL( std::fabs( podOutput->residualStandardDeviation_ ), 1.0E-3 );

    // Check parameter errors
    for( unsigned int i = 0; i < 4; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( i ) - truthParameters( i ) ), 1.0E-13 );
    }
    for( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( i + 4 ) - truthParameters( i + 4 ) ), 1.0E-17 );
    }
    for( unsigned int i = 0; i < 1; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( i + 7 ) - truthParameters( i + 7 ) ), 1.0E-16 );

    }
    std::cout<<"True error: "<<( podOutput->parameterEstimate_ - truthParameters ).transpose( )<<std::endl;
    std::cout<<"Formal error: "<<podOutput->getFormalErrorVector( ).transpose( )<<std::endl;
    std::cout<<"Error ratio: "<<( ( 1.0E-3 * podOutput->getFormalErrorVector( ).segment( 0, numberOfParameters ) ).cwiseQuotient(
                                      podOutput->parameterEstimate_ - truthParameters ) ).transpose( )<<std::endl;

//    input_output::writeMatrixToFile( podOutput->normalizedInformationMatrix_,
//                                     "rotationTestEstimationInformationMatrix.dat", 16,
//                                     input_output::getTudatRootPath( ) );
//    input_output::writeMatrixToFile( podOutput->informationMatrixTransformationDiagonal_,
//                                     "rotationTestEstimationInformationMatrixNormalization.dat", 16,
//                                     input_output::getTudatRootPath( ) );
//    input_output::writeMatrixToFile( podOutput->weightsMatrixDiagonal_,
//                                     "rotationTestEstimationWeightsDiagonal.dat", 16,
//                                     input_output::getTudatRootPath( ) );
//    input_output::writeMatrixToFile( podOutput->residuals_,
//                                     "rotationTestEstimationResiduals.dat", 16,
//                                     input_output::getTudatRootPath( ) );
//    input_output::writeMatrixToFile( podOutput->getCorrelationMatrix( ),
//                                     "rotationTestEstimationCorrelations.dat", 16,
//                                     input_output::getTudatRootPath( ) );
//    input_output::writeMatrixToFile( podOutput->getResidualHistoryMatrix( ),
//                                     "rotationTestResidualHistory.dat", 16,
//                                     input_output::getTudatRootPath( ) );
//    input_output::writeMatrixToFile( podOutput->getParameterHistoryMatrix( ),
//                                     "rotationTestParameterHistory.dat", 16,
//                                     input_output::getTudatRootPath( ) );
//    input_output::writeMatrixToFile( getConcatenatedMeasurementVector( podInput->getObservationsAndTimes( ) ),
//                                     "rotationTestObservationMeasurements.dat", 16,
//                                     input_output::getTudatRootPath( ) );
//    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
//                                         getConcatenatedTimeVector( podInput->getObservationsAndTimes( ) ) ),
//                                     "rotationTestObservationTimes.dat", 16,
//                                     input_output::getTudatRootPath( ) );
//    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
//                                         getConcatenatedGroundStationIndex( podInput->getObservationsAndTimes( ) ).first ),
//                                     "rotationTestObservationLinkEnds.dat", 16,
//                                     input_output::getTudatRootPath( ) );
//    input_output::writeMatrixToFile( getConcatenatedMeasurementVector( podInput->getObservationsAndTimes( ) ),
//                                     "rotationTestObservationMeasurements.dat", 16,
//                                     input_output::getTudatRootPath( ) );

//    input_output::writeMatrixToFile( podOutput->getFormalErrorVector( ),
//                                     "rotationTestObservationFormalEstimationError.dat", 16,
//                                     input_output::getTudatRootPath( ) );
}

//! Test if Phobos coupled translational-rotational dynamics is correctly estimation from lander tracking data
BOOST_AUTO_TEST_CASE( test_RotationalTranslationalDynamicsEstimationFromLanderData )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Retrieve list of body objects.
    SystemOfBodies bodies = getTestBodyMap( 9376.0E3, 0 );
    createGroundStation( bodies.at( "Phobos" ), "Lander", ( Eigen::Vector3d( ) << 0.1, 0.35, 0.0 ).finished( ), geodetic_position );

    // Define time range of test.
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 60.0 * 86400.0;

    // Define mean motion (equal to rotation rate).
    double phobosSemiMajorAxis = 9376.0E3;
    double meanMotion = std::sqrt( getBodyGravitationalParameter( "Mars" ) /
                                   std::pow( phobosSemiMajorAxis, 3.0 ) );

    // Define initial rotational state
    double initialAnglePerturbation = 1.0E-2;
    double initialRotationRatePerturbation = 1.0E-2;

    Eigen::Quaterniond nominalInitialRotation =
            Eigen::Quaterniond( Eigen::AngleAxisd( -initialAnglePerturbation, Eigen::Vector3d::UnitZ( ) ) );
    Eigen::VectorXd systemInitialState = Eigen::VectorXd::Zero( 7 );
    systemInitialState.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat( nominalInitialRotation );
    systemInitialState( 6 ) = meanMotion * ( 1.0 + initialRotationRatePerturbation );

    Eigen::Matrix3d phobosInertiaTensor = bodies.at( "Phobos" )->getBodyInertiaTensor( );

    // Create torque models
    SelectedTorqueMap torqueMap;
    torqueMap[ "Phobos" ][ "Mars" ].push_back( std::make_shared< TorqueSettings >( second_order_gravitational_torque ) );

    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "Phobos" );

    basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap(
                bodies, torqueMap, bodiesToIntegrate );

    // Define acceleration model settings.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > translationalBodiesToPropagate;
    std::vector< std::string > translationalCentralBodies;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfPhobos;
    accelerationsOfPhobos[ "Mars" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 2, 2, 2, 2 ) );
    accelerationMap[ "Phobos" ] = accelerationsOfPhobos;
    translationalBodiesToPropagate.push_back( "Phobos" );
    translationalCentralBodies.push_back( "Mars" );

    // Create acceleration models
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, translationalBodiesToPropagate, translationalCentralBodies );

    // Define integrator settings.
    double timeStep = 240.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettings< > >
            ( initialEphemerisTime, timeStep,
              rungeKuttaFehlberg78, timeStep, timeStep, 1.0, 1.0 );

    // Define propagator settings.
    std::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
            std::make_shared< RotationalStatePropagatorSettings< double > >
            ( torqueModelMap, bodiesToIntegrate, systemInitialState, std::make_shared< PropagationTimeTerminationSettings >(
                  finalEphemerisTime ) );

    Eigen::VectorXd initialTranslationalState =
            propagators::getInitialStatesOfBodies(
                translationalBodiesToPropagate, translationalCentralBodies, bodies, initialEphemerisTime );
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( translationalCentralBodies, accelerationModelMap, translationalBodiesToPropagate,
              initialTranslationalState, finalEphemerisTime );

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > >  propagatorSettingsList;
    propagatorSettingsList.push_back( translationalPropagatorSettings );
    propagatorSettingsList.push_back( rotationalPropagatorSettings );

    std::shared_ptr< PropagatorSettings< double > > propagatorSettings =
            std::make_shared< MultiTypePropagatorSettings< double > >(
                propagatorSettingsList,
                std::make_shared< PropagationTimeTerminationSettings >( finalEphemerisTime ) );

    // Define link ends and observables
    std::vector< LinkEnds > linkEndsList;
    LinkEnds currentLinkEnds;
    currentLinkEnds[ transmitter ] = std::make_pair( "Earth", "" );
    currentLinkEnds[ receiver ] = std::make_pair( "Phobos", "Lander" );
    linkEndsList.push_back( currentLinkEnds );

    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_range ].push_back( linkEndsList[ 0 ] );

    // Create parameters to estimate
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames = getInitialStateParameterSettings< double >( propagatorSettings, bodies );

    std::vector< std::pair< int, int > > blockIndices;
    blockIndices.push_back( std::make_pair( 2, 0 ) );
    blockIndices.push_back( std::make_pair( 2, 2 ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  blockIndices, "Phobos", spherical_harmonics_cosine_coefficient_block ) );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodies );
    printEstimatableParameterEntries( parametersToEstimate );

    // Create observation settings
    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;

        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            observationSettingsList.push_back(
                        std::make_shared< ObservationModelSettings >( currentObservable, currentLinkEndsList.at( i ) ) );
        }
    }

    // Create orbit determination object
    OrbitDeterminationManager< double, double > orbitDeterminationManager =
            OrbitDeterminationManager< double, double >(
                bodies, parametersToEstimate, observationSettingsList,
                integratorSettings, propagatorSettings );

    // Deifne observation times
    std::vector< double > observationTimes;
    double currentTime = initialEphemerisTime + 1800.0;
    double observationTimeStep = 1200.0;
    while( currentTime < finalEphemerisTime - 1800.0 )
    {
        observationTimes.push_back( currentTime );
        currentTime += observationTimeStep;
    }


    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            measurementSimulationInput.push_back(
                    std::make_shared< TabulatedObservationSimulationSettings< > >(
                        currentObservable, currentLinkEndsList.at( i ), observationTimes, receiver ) );
        }
    }

    // Simulate observations
    std::shared_ptr< ObservationCollection< > > observationsAndTimes = simulateObservations< double, double >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    // Perturb parameter estimate
    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< double >( );
    int numberOfParameters = initialParameterEstimate.rows( );

    Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    initialParameterEstimate( 7 ) -= 1.0E-5;
    initialParameterEstimate( 8 ) += 1.0E-5;
    initialParameterEstimate( 9 ) += 1.0E-5;

    initialParameterEstimate.segment( 6, 4 ).normalize( );
    initialParameterEstimate( 10 ) += 1.0E-7;

    initialParameterEstimate( 0 ) += 1.0;
    initialParameterEstimate( 1 ) += 1.0;
    initialParameterEstimate( 2 ) += 1.0;
    int parameterSize = initialParameterEstimate.rows( );

    // Define estimation input
    std::shared_ptr< PodInput< double, double  > > podInput =
            std::make_shared< PodInput< double, double > >(
                observationsAndTimes, parameterSize,
                Eigen::MatrixXd::Zero( parameterSize, parameterSize ), initialParameterEstimate - truthParameters );
    podInput->defineEstimationSettings( true, false, true, true, true, true );

    // Perform estimation
    std::shared_ptr< PodOutput< double > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, std::make_shared< EstimationConvergenceChecker >( 6 ) );

    // Check residual size (sub-mm over >1 AU)
    BOOST_CHECK_SMALL( std::fabs( podOutput->residualStandardDeviation_ ), 1.0E-3 );

    // Check parameter errors
    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 0 ) - truthParameters( 0 ) ), 1.0E-4 );
    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 1 ) - truthParameters( 1 ) ), 1.0E-4 );
    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 2 ) - truthParameters( 2 ) ), 1.0E-2 );

    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 3 ) - truthParameters( 3 ) ), 2.0E-8 );
    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 4 ) - truthParameters( 4 ) ), 2.0E-8 );
    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 5 ) - truthParameters( 5 ) ), 2.0E-6 );

    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 6 ) - truthParameters( 6 ) ), 1.0E-12 );
    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 7 ) - truthParameters( 7 ) ), 1.0E-9 );
    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 8 ) - truthParameters( 8 ) ), 1.0E-9 );
    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 9 ) - truthParameters( 9 ) ), 1.0E-9 );

    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 10 ) - truthParameters( 10 ) ), 5.0E-13 );
    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 11 ) - truthParameters( 11 ) ), 5.0E-13 );
    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 12 ) - truthParameters( 12 ) ), 5.0E-13 );

    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 13 ) - truthParameters( 13 ) ), 1.0E-10 );
    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 14 ) - truthParameters( 14 ) ), 1.0E-11 );
}

BOOST_AUTO_TEST_SUITE_END( )

}

}




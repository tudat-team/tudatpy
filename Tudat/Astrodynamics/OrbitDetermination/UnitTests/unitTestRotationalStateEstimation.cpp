/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/Ephemerides/keplerEphemeris.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"
#include "Tudat/SimulationSetup/tudatEstimationHeader.h"

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

NamedBodyMap getTestBodyMap( const double phobosSemiMajorAxis,
                             const bool useSymmetricEquator = 0 )
{
    NamedBodyMap bodyMap;
    bodyMap[ "Mars" ] = boost::make_shared< Body >( );
    bodyMap[ "Mars" ]->setEphemeris( createBodyEphemeris(
                                         getDefaultEphemerisSettings( "Mars" ), "Mars" ) );
    boost::shared_ptr< SphericalHarmonicsGravityFieldSettings > marsGravityFieldSettings =
            boost::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >(
                getDefaultGravityFieldSettings( "Mars", TUDAT_NAN, TUDAT_NAN ) );
    marsGravityFieldSettings->resetGravitationalParameter( spice_interface::getBodyGravitationalParameter( "Mars" ) );
    bodyMap[ "Mars" ]->setGravityFieldModel(
                createGravityFieldModel( marsGravityFieldSettings, "Mars", bodyMap ) );
//                boost::make_shared< gravitation::GravityFieldModel >(
//                    spice_interface::getBodyGravitationalParameter( "Mars" ) ) );

    bodyMap[ "Mars" ]->setRotationalEphemeris(
                createRotationModel( getDefaultRotationModelSettings(
                                         "Mars", TUDAT_NAN, TUDAT_NAN ), "Mars" ) );

    bodyMap[ "Earth" ] = boost::make_shared< Body >( );
    bodyMap[ "Earth" ]->setEphemeris( createBodyEphemeris(
                                          getDefaultEphemerisSettings( "Earth" ), "Earth" ) );

    bodyMap[ "Phobos" ] = boost::make_shared< Body >( );

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
    bodyMap[ "Phobos" ]->setBodyInertiaTensor( phobosInertiaTensor );
    bodyMap[ "Phobos" ]->setShapeModel(
                boost::make_shared< SphericalBodyShapeModel >( 15.0E3 ) );

    double phobosGravitationalParameter = phobosMass * physical_constants::GRAVITATIONAL_CONSTANT;

    Eigen::MatrixXd phobosCosineGravityFieldCoefficients = Eigen::Matrix3d::Zero( ),
            phobosSineGravityFieldCoefficients = Eigen::Matrix3d::Zero( );
    double phobosScaledMeanMomentOfInertia;
    gravitation::getDegreeTwoSphericalHarmonicCoefficients(
                phobosInertiaTensor, phobosGravitationalParameter, phobosReferenceRadius, true,
                phobosCosineGravityFieldCoefficients, phobosSineGravityFieldCoefficients, phobosScaledMeanMomentOfInertia );

    bodyMap[ "Phobos" ]->setGravityFieldModel(
                boost::make_shared< gravitation::SphericalHarmonicsGravityField >(
                    phobosGravitationalParameter, phobosReferenceRadius, phobosCosineGravityFieldCoefficients,
                    phobosSineGravityFieldCoefficients, "Phobos_Fixed",
                     boost::bind( &Body::setBodyInertiaTensorFromGravityFieldAndExistingMeanMoment, bodyMap.at( "Phobos" ), true ) ) );



    Eigen::Quaterniond noRotationQuaternion = Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );
    Eigen::Matrix< double, 7, 1 > unitRotationState = Eigen::Matrix< double, 7, 1 >::Zero( );
    unitRotationState( 0 ) = noRotationQuaternion.w( );
    unitRotationState( 1 ) = noRotationQuaternion.x( );
    unitRotationState( 2 ) = noRotationQuaternion.y( );
    unitRotationState( 3 ) = noRotationQuaternion.z( );

    std::map< double, Eigen::Matrix< double, 7, 1 > > dummyRotationMap;
    dummyRotationMap[ -1.0E100 ] = unitRotationState;
    dummyRotationMap[ 1.0E100 ] = unitRotationState;

    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 7, 1 > > > dummyInterpolator =
            boost::make_shared< interpolators::LinearInterpolator< double, Eigen::Matrix< double, 7, 1 > > >( dummyRotationMap );
    bodyMap[ "Phobos" ]->setRotationalEphemeris( boost::make_shared< TabulatedRotationalEphemeris< double, double > >(
                                                     dummyInterpolator, "ECLIPJ2000", "Phobos_Fixed" ) );

    Eigen::Vector6d phobosKeplerElements = Eigen::Vector6d::Zero( );
    phobosKeplerElements( 0 ) = phobosSemiMajorAxis;

    bodyMap[ "Phobos" ]->setEphemeris(
                ephemerides::getTabulatedEphemeris( boost::make_shared< ephemerides::KeplerEphemeris >(
                                                        phobosKeplerElements, 0.0, spice_interface::getBodyGravitationalParameter( "Mars" ),
                                                        "Mars", "ECLIPJ2000" ), -3600.0, 40.0 * 86400.0 + 3600.0, 120.0,
                                                    boost::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ), "Mars" ) );

    return bodyMap;
}

BOOST_AUTO_TEST_CASE( test_RotationalDynamicsEstimationFromLanderData )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Retrieve list of body objects.
    NamedBodyMap bodyMap = getTestBodyMap( 9376.0E3, 0 );
    createGroundStation( bodyMap.at( "Phobos" ), "Lander", ( Eigen::Vector3d( ) << 0.1, 0.35, 0.0 ).finished( ), geodetic_position );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Define time range of test.
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 40.0 * 86400.0;

    // Set torques between bodies that are to be taken into account.
    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "Phobos" );

    // Define mean motion (equal to rotation rate).
    double phobosSemiMajorAxis = 9376.0E3;
    double meanMotion = std::sqrt( getBodyGravitationalParameter( "Mars" ) /
                                   std::pow( phobosSemiMajorAxis, 3.0 ) );

    // Define initial rotational state
    double initialAnglePerturbation = 1.0E-6;
    double initialRotationRatePerturbation = 1.0E-6;

    Eigen::Quaterniond nominalInitialRotation =
            Eigen::Quaterniond( Eigen::AngleAxisd( -initialAnglePerturbation, Eigen::Vector3d::UnitZ( ) ) );
    Eigen::VectorXd systemInitialState = Eigen::VectorXd::Zero( 7 );
    systemInitialState.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat( nominalInitialRotation );
    systemInitialState( 6 ) = meanMotion * ( 1.0 + initialRotationRatePerturbation );

    Eigen::Matrix3d phobosInertiaTensor = bodyMap.at( "Phobos" )->getBodyInertiaTensor( );

    // Create torque models
    SelectedTorqueMap torqueMap;
    torqueMap[ "Phobos" ][ "Mars" ].push_back( boost::make_shared< TorqueSettings >( second_order_gravitational_torque ) );

    // Create torque models
    basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap(
                bodyMap, torqueMap, bodiesToIntegrate );

    // Define integrator settings.
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, initialEphemerisTime, 5.0 );

    // Define propagator settings.
    boost::shared_ptr< RotationalStatePropagatorSettings< double > > propagatorSettings =
            boost::make_shared< RotationalStatePropagatorSettings< double > >
            ( torqueModelMap, bodiesToIntegrate, systemInitialState, boost::make_shared< PropagationTimeTerminationSettings >(
                  finalEphemerisTime ) );

    std::vector< LinkEnds > linkEndsList;
    LinkEnds currentLinkEnds;
    currentLinkEnds[ transmitter ] = std::make_pair( "Earth", "" );
    currentLinkEnds[ receiver ] = std::make_pair( "Phobos", "Lander" );
    linkEndsList.push_back( currentLinkEnds );

    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_range ].push_back( linkEndsList[ 0 ] );

    std::vector< boost::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back(
                boost::make_shared< InitialRotationalStateEstimatableParameterSettings< double > >(
                    "Phobos", systemInitialState ) );
    parameterNames.push_back( boost::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  2, 0, 2, 2, "Phobos", spherical_harmonics_cosine_coefficient_block ) );

    // Create parameters
    boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodyMap );
    printEstimatableParameterEntries( parametersToEstimate );

    observation_models::ObservationSettingsMap observationSettingsMap;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;

        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            observationSettingsMap.insert( std::make_pair( currentLinkEndsList.at( i ),
                                                           boost::make_shared< ObservationSettings >( currentObservable ) ) );
        }
    }

//    std::cout<<"Pre inertia tensor "<<std::setprecision( 16 )<<bodyMap.at( "Phobos" )->getBodyInertiaTensor( )<<std::endl;
//    std::cout<<std::setprecision( 6 )<<std::endl;

    OrbitDeterminationManager< double, double > orbitDeterminationManager =
            OrbitDeterminationManager< double, double >(
                bodyMap, parametersToEstimate, observationSettingsMap,
                integratorSettings, propagatorSettings );

    std::vector< double > observationTimes;
    double currentTime = initialEphemerisTime + 1800.0;
    double observationTimeStep = 120.0;
    while( currentTime < finalEphemerisTime - 1800.0 )
    {
        observationTimes.push_back( currentTime );
        currentTime += observationTimeStep;
    }

    std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< double >, LinkEndType > > > measurementSimulationInput;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            measurementSimulationInput[ currentObservable ][ currentLinkEndsList.at( i ) ] =
                    std::make_pair( observationTimes, receiver );
        }
    }

    typedef Eigen::Matrix< double, Eigen::Dynamic, 1 > ObservationVectorType;
    typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< double >, LinkEndType > > > SingleObservablePodInputType;
    typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;

    // Simulate observations
    PodInputDataType observationsAndTimes = simulateObservations< double, double >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ) );

    // Perturb parameter estimate
    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< double >( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    int numberOfParameters = initialParameterEstimate.rows( );
//    initialParameterEstimate( 1 ) += 1.0E-5;
    initialParameterEstimate( 2 ) -= 1.0E-5;

    initialParameterEstimate.segment( 0, 4 ).normalize( );
    initialParameterEstimate( 4 ) += 1.0E-7;

    initialParameterEstimate( 8 ) += 1.0E-3;
    initialParameterEstimate( 9 ) += 1.0E-3;

    // Define estimation input
    boost::shared_ptr< PodInput< double, double  > > podInput =
            boost::make_shared< PodInput< double, double > >(
                observationsAndTimes, numberOfParameters,
                Eigen::MatrixXd::Zero( numberOfParameters, numberOfParameters ), initialParameterEstimate - truthParameters );

    // Perform estimation
    boost::shared_ptr< PodOutput< double > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, boost::make_shared< EstimationConvergenceChecker >( 3 ) );

//    std::cout<<podOutput->parameterEstimate_.transpose( )<<std::endl;
    for( unsigned int i = 0; i < 4; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( i ) - truthParameters( i ) ), 1.0E-10 );
    }

    for( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( i + 4 ) - truthParameters( i + 4 ) ), 1.0E-14 );
    }
    std::cout<<( podOutput->parameterEstimate_ - truthParameters ).transpose( )<<std::endl;

//    std::cout<<"Post inertia tensor "<<std::setprecision( 16 )<<bodyMap.at( "Phobos" )->getBodyInertiaTensor( )<<std::endl;

    input_output::writeMatrixToFile( podOutput->normalizedInformationMatrix_,
                                     "rotationInformationMatrix.dat", 16 );
    input_output::writeMatrixToFile( podOutput->getCorrelationMatrix( ),
                                     "rotationCorrelations.dat", 16 );
    input_output::writeMatrixToFile( podOutput->inverseNormalizedCovarianceMatrix_,
                                     "rotationInverseNormalizedCovariance.dat", 16 );
    input_output::writeMatrixToFile( podOutput->getFormalErrorVector( ),
                                     "rotationFormalEstimationError.dat", 16 );
    input_output::writeMatrixToFile( truthParameters,
                                     "rotationTruthParameters.dat", 16 );
    input_output::writeMatrixToFile( podOutput->residuals_,
                                     "rotationResiduals.dat", 16 );
    input_output::writeMatrixToFile( podOutput->informationMatrixTransformationDiagonal_,
                                     "rotationParameterNormalization.dat", 16 );
}

BOOST_AUTO_TEST_CASE( test_RotationalTranslationalDynamicsEstimationFromLanderData )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Retrieve list of body objects.
    NamedBodyMap bodyMap = getTestBodyMap( 9376.0E3, 0 );
    createGroundStation( bodyMap.at( "Phobos" ), "Lander", ( Eigen::Vector3d( ) << 0.1, 0.35, 0.0 ).finished( ), geodetic_position );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Define time range of test.
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 60.0 * 86400.0;

    // Set torques between bodies that are to be taken into account.
    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "Phobos" );

    // Define mean motion (equal to rotation rate).
    double phobosSemiMajorAxis = 9376.0E3;
    double meanMotion = std::sqrt( getBodyGravitationalParameter( "Mars" ) /
                                   std::pow( phobosSemiMajorAxis, 3.0 ) );

    // Define initial rotational state
    double initialAnglePerturbation = 1.0E-6;
    double initialRotationRatePerturbation = 1.0E-6;

    Eigen::Quaterniond nominalInitialRotation =
            Eigen::Quaterniond( Eigen::AngleAxisd( -initialAnglePerturbation, Eigen::Vector3d::UnitZ( ) ) );
    Eigen::VectorXd systemInitialState = Eigen::VectorXd::Zero( 7 );
    systemInitialState.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat( nominalInitialRotation );
    systemInitialState( 6 ) = meanMotion * ( 1.0 + initialRotationRatePerturbation );

    Eigen::Matrix3d phobosInertiaTensor = bodyMap.at( "Phobos" )->getBodyInertiaTensor( );

    // Create torque models
    SelectedTorqueMap torqueMap;
    torqueMap[ "Phobos" ][ "Mars" ].push_back( boost::make_shared< TorqueSettings >( second_order_gravitational_torque ) );

    // Create torque models
    basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap(
                bodyMap, torqueMap, bodiesToIntegrate );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > translationalBodiesToPropagate;
    std::vector< std::string > translationalCentralBodies;

    // Define acceleration model settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfPhobos;
    accelerationsOfPhobos[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    //accelerationsOfPhobos[ "Mars" ].push_back( boost::make_shared< MutualSphericalHarmonicAccelerationSettings >( 7, 7, 2, 2 ) );

    accelerationMap[ "Phobos" ] = accelerationsOfPhobos;

    translationalBodiesToPropagate.push_back( "Phobos" );
    translationalCentralBodies.push_back( "Mars" );

    // Create acceleration models
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, translationalBodiesToPropagate, translationalCentralBodies );

    // Define integrator settings.
    double timeStep = 180.0;
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
            ( rungeKuttaVariableStepSize, initialEphemerisTime, timeStep,
              RungeKuttaCoefficients::rungeKuttaFehlberg78, timeStep, timeStep, 1.0, 1.0 );


    // Define propagator settings.
    boost::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
            boost::make_shared< RotationalStatePropagatorSettings< double > >
            ( torqueModelMap, bodiesToIntegrate, systemInitialState, boost::make_shared< PropagationTimeTerminationSettings >(
                  finalEphemerisTime ) );

    Eigen::VectorXd initialTranslationalState =
            propagators::getInitialStatesOfBodies(
                translationalBodiesToPropagate, translationalCentralBodies, bodyMap, initialEphemerisTime );
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( translationalCentralBodies, accelerationModelMap, translationalBodiesToPropagate,
              initialTranslationalState, finalEphemerisTime );

    std::vector< boost::shared_ptr< SingleArcPropagatorSettings< double > > >  propagatorSettingsList;
    propagatorSettingsList.push_back( translationalPropagatorSettings );
    propagatorSettingsList.push_back( rotationalPropagatorSettings );

    boost::shared_ptr< PropagatorSettings< double > > propagatorSettings =
            boost::make_shared< MultiTypePropagatorSettings< double > >(
                propagatorSettingsList,
                boost::make_shared< PropagationTimeTerminationSettings >( finalEphemerisTime ) );


    std::vector< LinkEnds > linkEndsList;
    LinkEnds currentLinkEnds;
    currentLinkEnds[ transmitter ] = std::make_pair( "Earth", "" );
    currentLinkEnds[ receiver ] = std::make_pair( "Phobos", "Lander" );
    linkEndsList.push_back( currentLinkEnds );

    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_range ].push_back( linkEndsList[ 0 ] );

    std::vector< boost::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back(
                boost::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                    "Phobos", initialTranslationalState, "Mars" ) );
    parameterNames.push_back(
                boost::make_shared< InitialRotationalStateEstimatableParameterSettings< double > >(
                    "Phobos", systemInitialState ) );


    // Create parameters
    boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodyMap );
    printEstimatableParameterEntries( parametersToEstimate );

    observation_models::ObservationSettingsMap observationSettingsMap;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;

        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            observationSettingsMap.insert( std::make_pair( currentLinkEndsList.at( i ),
                                                           boost::make_shared< ObservationSettings >( currentObservable ) ) );
        }
    }

    OrbitDeterminationManager< double, double > orbitDeterminationManager =
            OrbitDeterminationManager< double, double >(
                bodyMap, parametersToEstimate, observationSettingsMap,
                integratorSettings, propagatorSettings );

    std::vector< double > observationTimes;
    double currentTime = initialEphemerisTime + 1800.0;
    double observationTimeStep = 1200.0;
    while( currentTime < finalEphemerisTime - 1800.0 )
    {
        observationTimes.push_back( currentTime );
        currentTime += observationTimeStep;
    }

    std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< double >, LinkEndType > > > measurementSimulationInput;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            measurementSimulationInput[ currentObservable ][ currentLinkEndsList.at( i ) ] =
                    std::make_pair( observationTimes, receiver );
        }
    }

    typedef Eigen::Matrix< double, Eigen::Dynamic, 1 > ObservationVectorType;
    typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< double >, LinkEndType > > > SingleObservablePodInputType;
    typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;

    // Simulate observations
    PodInputDataType observationsAndTimes = simulateObservations< double, double >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ) );

    // Perturb parameter estimate
    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< double >( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    //    initialParameterEstimate( 1 ) += 1.0E-5;
    initialParameterEstimate( 8 ) -= 1.0E-5;

    initialParameterEstimate.segment( 6, 4 ).normalize( );
    initialParameterEstimate( 10 ) += 1.0E-7;

    initialParameterEstimate( 0 ) += 1.0;
    initialParameterEstimate( 1 ) += 1.0;
    initialParameterEstimate( 2 ) += 1.0;

    std::cout<<"Truth: "<<( truthParameters ).transpose( )<<std::endl;
    std::cout<<"New: "<<( initialParameterEstimate ).transpose( )<<std::endl;

    // Define estimation input
    boost::shared_ptr< PodInput< double, double  > > podInput =
            boost::make_shared< PodInput< double, double > >(
                observationsAndTimes, 13,
                Eigen::MatrixXd::Zero( 13, 13 ), initialParameterEstimate - truthParameters );
    podInput->defineEstimationSettings( true, false );

    // Perform estimation
    boost::shared_ptr< PodOutput< double > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, boost::make_shared< EstimationConvergenceChecker >( 4 ) );

    //    std::cout<<podOutput->parameterEstimate_.transpose( )<<std::endl;
    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 0 ) - truthParameters( 0 ) ), 1.0E-4 );
    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 1 ) - truthParameters( 1 ) ), 1.0E-2 );
    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 2 ) - truthParameters( 2 ) ), 1.0E-2 );

    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 3 ) - truthParameters( 3 ) ), 1.0E-6 );
    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 4 ) - truthParameters( 4 ) ), 1.0E-8 );
    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 5 ) - truthParameters( 5 ) ), 1.0E-4 );

    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 6 ) - truthParameters( 6 ) ), 1.0E-14 );
    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 7 ) - truthParameters( 7 ) ), 1.0E-8 );
    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 8 ) - truthParameters( 8 ) ), 1.0E-8 );
    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 9 ) - truthParameters( 9 ) ), 1.0E-8 );

    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 10 ) - truthParameters( 10 ) ), 1.0E-12 );
    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 11 ) - truthParameters( 11 ) ), 1.0E-12 );
    BOOST_CHECK_SMALL( std::fabs( podOutput->parameterEstimate_( 12 ) - truthParameters( 12 ) ), 1.0E-12 );

    std::cout<<( podOutput->parameterEstimate_ - truthParameters ).transpose( )<<std::endl;

    input_output::writeMatrixToFile( podOutput->normalizedInformationMatrix_,
                                     "rotationInformationMatrix.dat", 16 );
    input_output::writeMatrixToFile( podOutput->getCorrelationMatrix( ),
                                     "rotationCorrelations.dat", 16 );
    input_output::writeMatrixToFile( podOutput->inverseNormalizedCovarianceMatrix_,
                                     "rotationInverseNormalizedCovariance.dat", 16 );
    input_output::writeMatrixToFile( podOutput->getFormalErrorVector( ),
                                     "rotationFormalEstimationError.dat", 16 );
    input_output::writeMatrixToFile( truthParameters,
                                     "rotationTruthParameters.dat", 16 );
    input_output::writeMatrixToFile( podOutput->residuals_,
                                     "rotationResiduals.dat", 16 );
    input_output::writeMatrixToFile( podOutput->informationMatrixTransformationDiagonal_,
                                     "rotationParameterNormalization.dat", 16 );

}

BOOST_AUTO_TEST_SUITE_END( )

}

}




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

#include "tudat/simulation/estimation_setup/orbitDeterminationTestCases.h"
#include "tudat/simulation/estimation_setup/podProcessing.h"
#include "tudat/astro/propagators/propagateCovariance.h"

namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_covariance_propagation )

BOOST_AUTO_TEST_CASE( test_EstimationInputAndOutput )
{
    const double startTime = 1.0E7;

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    for( int test = 0; test < 2; test++ )
    {
        // Define bodies in simulation
        std::vector< std::string > bodyNames;
        bodyNames.push_back( "Earth" );
        bodyNames.push_back( "Sun" );
        bodyNames.push_back( "Moon" );
        bodyNames.push_back( "Mars" );

        // Specify initial time
        double initialEphemerisTime = startTime;
        double finalEphemerisTime = initialEphemerisTime + 4.0 * 3600.0;

        // Create bodies needed in simulation
        BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
        bodySettings.at( "Earth" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
            "ECLIPJ2000", "IAU_Earth",
            spice_interface::computeRotationQuaternionBetweenFrames(
                "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
            initialEphemerisTime, 2.0 * mathematical_constants::PI /
                                  ( physical_constants::JULIAN_DAY ) );

        SystemOfBodies bodies = createSystemOfBodies( bodySettings );
        bodies.createEmptyBody( "Vehicle" );
        bodies.at( "Vehicle" )->setConstantBodyMass( 400.0 );

        // Create aerodynamic coefficient interface settings.
        double referenceArea = 4.0;
        double aerodynamicCoefficient = 1.2;
        std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >(
                referenceArea, aerodynamicCoefficient * ( Eigen::Vector3d( ) << 1.2, -0.1, -0.4 ).finished( ),
                negative_aerodynamic_frame_coefficients );

        // Create and set aerodynamic coefficnumberOfDaysOfDataients object
        bodies.at( "Vehicle" )->setAerodynamicCoefficientInterface(
            createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle", bodies ) );

        // Create radiation pressure settings
        double referenceAreaRadiation = 4.0;
        double radiationPressureCoefficient = 1.2;
        std::vector< std::string > occultingBodies;
        occultingBodies.push_back( "Earth" );
        std::shared_ptr< RadiationPressureInterfaceSettings > asterixRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

        // Create and set radiation pressure settings
        bodies.at( "Vehicle" )->setRadiationPressureInterface(
            "Sun", createRadiationPressureInterface(
                asterixRadiationPressureSettings, "Vehicle", bodies ) );

    //    bodies.at( "Vehicle" )->setEphemeris( std::make_shared< TabulatedCartesianEphemeris< > >(
    //        std::shared_ptr< interpolators::OneDimensionalInterpolator
    //            < double, Eigen::Vector6d > >( ), "Earth", "ECLIPJ2000" ) );

        // Set accelerations on Vehicle that are to be taken into account.
        SelectedAccelerationMap accelerationMap;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
        accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );
        accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
            basic_astrodynamics::point_mass_gravity ) );
        accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
            basic_astrodynamics::point_mass_gravity ) );
        accelerationsOfVehicle[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
            basic_astrodynamics::point_mass_gravity ) );
        accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
            basic_astrodynamics::cannon_ball_radiation_pressure ) );
        accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
            basic_astrodynamics::aerodynamic ) );
        accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

        // Set bodies for which initial state is to be estimated and integrated.
        std::vector< std::string > bodiesToIntegrate;
        std::vector< std::string > centralBodies;
        bodiesToIntegrate.push_back( "Vehicle" );
        centralBodies.push_back( "Earth" );

        // Create acceleration models
        AccelerationMap accelerationModelMap = createAccelerationModelsMap(
            bodies, accelerationMap, bodiesToIntegrate, centralBodies );

        // Set Keplerian elements for Asterix.
        Eigen::Vector6d asterixInitialStateInKeplerianElements;
        asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7200.0E3;
        asterixInitialStateInKeplerianElements( eccentricityIndex ) = 1.0E-5;
        asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
        asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = unit_conversions::convertDegreesToRadians( 235.7 );
        asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 23.4 );
        asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

        double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

        // Set (perturbed) initial state.
        Eigen::Matrix< double, 6, 1 > systemInitialState = convertKeplerianToCartesianElements(
            asterixInitialStateInKeplerianElements, earthGravitationalParameter );

        // Create propagator settings
        std::shared_ptr< IntegratorSettings< double > > integratorSettings =
            rungeKuttaFixedStepSettings( 400.0, CoefficientSets::rungeKuttaFehlberg78 );
        std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double, double > >
                ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, initialEphemerisTime,
                  integratorSettings, propagationTimeTerminationSettings( finalEphemerisTime ), cowell );


        // Define parameters.
        LinkEnds linkEnds;
        linkEnds[ observed_body ] = LinkEndId( "Vehicle", "" );

        std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
        linkEndsPerObservable[ position_observable ].push_back( linkEnds );

        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
        parameterNames.push_back(
            std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                "Vehicle", systemInitialState, "Earth" ) );

        if( test > 0 )
        {
            parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", radiation_pressure_coefficient ) );
            parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", constant_drag_coefficient ) );
        }

        // Create parameters
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, bodies );

        printEstimatableParameterEntries( parametersToEstimate );

        std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
            position_observable, linkEnds ) );

        // Create orbit determination object.
        OrbitDeterminationManager< double, double > orbitDeterminationManager =
            OrbitDeterminationManager< double, double >(
                bodies, parametersToEstimate, observationSettingsList, propagatorSettings );

        std::vector< double > baseTimeList;
        double observationTimeStart = initialEphemerisTime + 100.0;
        double observationInterval = 20.0;
        double currentObservationTime = observationTimeStart;
        while( currentObservationTime < finalEphemerisTime - 100.0 )
        {
            baseTimeList.push_back( currentObservationTime );
            currentObservationTime += observationInterval;
        }

        std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput =
            getObservationSimulationSettings< double >(
                linkEndsPerObservable, baseTimeList, observed_body );

        // Simulate observations
        std::shared_ptr< ObservationCollection< double, double > > simulatedObservations =
            simulateObservations< double, double >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );



        // Define estimation input
        std::shared_ptr< CovarianceAnalysisInput< double, double  > > estimationInput =
            std::make_shared< CovarianceAnalysisInput< double, double > >(
                simulatedObservations );

        std::map< observation_models::ObservableType, double > weightPerObservable;
        weightPerObservable[ position_observable ] = 1.0;
        estimationInput->setConstantPerObservableWeightsMatrix( weightPerObservable );

        // Perform estimation
        std::shared_ptr< CovarianceAnalysisOutput< double > > estimationOutput = orbitDeterminationManager.computeCovariance(
            estimationInput );
        auto stateHistory = std::dynamic_pointer_cast< SingleArcVariationalEquationsSolver< > >( orbitDeterminationManager.getVariationalEquationsSolver( )
        )->getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
        auto stateTransitionMatrixHistory =
            std::dynamic_pointer_cast< SingleArcVariationalEquationsSolver< > >( orbitDeterminationManager.getVariationalEquationsSolver( )
            )->getStateTransitionMatrixSolution( );
        auto sensitivityMatrixHistory =
            std::dynamic_pointer_cast< SingleArcVariationalEquationsSolver< > >( orbitDeterminationManager.getVariationalEquationsSolver( )
            )->getSensitivityMatrixSolution( );

        std::map< double, Eigen::MatrixXd > propagatedCovariance;
        propagateCovariance(
            propagatedCovariance,
            estimationOutput->getUnnormalizedCovarianceMatrix( ),
            orbitDeterminationManager.getStateTransitionAndSensitivityMatrixInterface( ),
            utilities::createVectorFromMapKeys( stateTransitionMatrixHistory ) );

        std::map< double, Eigen::MatrixXd > rswPropagatedCovariance = convertCovarianceHistoryToFrame(
            propagatedCovariance,
            stateHistory,
            reference_frames::global_reference_frame,
            reference_frames::rsw_reference_frame );
        std::map< double, Eigen::MatrixXd > tnwPropagatedCovariance = convertCovarianceHistoryToFrame(
            propagatedCovariance,
            stateHistory,
            reference_frames::global_reference_frame,
            reference_frames::tnw_reference_frame );

        for( auto it : propagatedCovariance )
        {
            Eigen::MatrixXd currentStateTransition = stateTransitionMatrixHistory.at( it.first );
            if( test == 0 )
            {
                Eigen::MatrixXd manualCovariance = currentStateTransition * estimationOutput->getUnnormalizedCovarianceMatrix( ) * currentStateTransition.transpose( );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( it.second,  manualCovariance, 1.0E-10 );

                Eigen::Matrix3d rotationToRsw = reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrix( stateHistory.at( it.first ) );
                Eigen::Matrix6d rotationToRswFull = Eigen::Matrix6d::Zero( );
                rotationToRswFull.block( 0, 0, 3, 3 ) = rotationToRsw;
                rotationToRswFull.block( 3, 3, 3, 3 ) = rotationToRsw;
                Eigen::MatrixXd manualRswCovariance = rotationToRswFull * it.second * rotationToRswFull.transpose( );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( rswPropagatedCovariance.at( it.first ),  manualRswCovariance, 1.0E-9 );

                Eigen::Matrix3d rotationToTnw = reference_frames::getInertialToTnwRotation( stateHistory.at( it.first ) );
                Eigen::Matrix6d rotationToTnwFull = Eigen::Matrix6d::Zero( );
                rotationToTnwFull.block( 0, 0, 3, 3 ) = rotationToTnw;
                rotationToTnwFull.block( 3, 3, 3, 3 ) = rotationToTnw;
                Eigen::MatrixXd manualTnwCovariance = rotationToTnwFull * it.second * rotationToTnwFull.transpose( );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( tnwPropagatedCovariance.at( it.first ),  manualTnwCovariance, 1.0E-9 );

                Eigen::VectorXd rswFormalError = rswPropagatedCovariance.at( it.first ).diagonal( ).cwiseSqrt( );
                Eigen::VectorXd tnwFormalError = tnwPropagatedCovariance.at( it.first ).diagonal( ).cwiseSqrt( );

                // Weak check on radial and along-track errors.
                BOOST_CHECK_SMALL( rswFormalError( 0 ) / rswFormalError( 1 ), 0.6 );
                BOOST_CHECK_SMALL( rswFormalError( 4 ) / rswFormalError( 3 ), 0.6 );

                // Check if in-plane RSW and TNW are similar
                BOOST_CHECK_CLOSE_FRACTION( rswFormalError( 0 ), tnwFormalError( 1 ), 1.0E-3 );
                BOOST_CHECK_CLOSE_FRACTION( rswFormalError( 1 ), tnwFormalError( 0 ), 1.0E-3 );

                BOOST_CHECK_CLOSE_FRACTION( rswFormalError( 3 ), tnwFormalError( 4 ), 1.0E-3 );
                BOOST_CHECK_CLOSE_FRACTION( rswFormalError( 4 ), tnwFormalError( 3 ), 1.0E-3 );

                // Check if cross-plane RSW and TNW are equal
                BOOST_CHECK_CLOSE_FRACTION( rswFormalError( 2 ), tnwFormalError( 2 ), std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_CLOSE_FRACTION( rswFormalError( 5 ), tnwFormalError( 5 ), std::numeric_limits< double >::epsilon( ) );
            }
            else
            {
                Eigen::MatrixXd currentSensitivity = sensitivityMatrixHistory.at( it.first );
                Eigen::MatrixXd fullStateTransition = Eigen::MatrixXd::Zero( 6, 8 );
                fullStateTransition.block( 0, 0, 6, 6 ) = currentStateTransition;
                fullStateTransition.block( 0, 6, 6, 2 ) = currentSensitivity;
                Eigen::MatrixXd manualCovariance = fullStateTransition * estimationOutput->getUnnormalizedCovarianceMatrix( ) * fullStateTransition.transpose( );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( it.second,  manualCovariance, 1.0E-12 );

            }
        }
//        std::cout <<"formal error: "<< ( estimationOutput->getFormalErrorVector( ) ).transpose( ) << std::endl;
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}



